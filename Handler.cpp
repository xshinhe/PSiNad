#include "Handler.h"

#include "core/DataSet.h"
#include "kernels/Kernel_Record.h"
#include "models/ModelFactory.h"
#include "solvers/SolverFactory.h"

//
#include "mpi_utils.h"


namespace PROJECT_NS {


Handler::Handler(hdlr_t itype, const std::string& solver_name, const std::string& model_name) {
    type   = itype;
    model  = ModelFactory(model_name);
    solver = SolverFactory(solver_name, model);
};


int Handler::run(Param* P) {
    switch (type) {
        case hdlr_t::single: {
            run_single(P);
            break;
        }
        case hdlr_t::multiple: {
            run_multiple(P);
            break;
        }
    }
    return 0;
}

int Handler::run_single(Param* P) {
    std::cout << P->repr();
    DataSet DS;
    auto begin = std::chrono::steady_clock::now();
    {
        solver->read_param(P);   // parameter parser and secondary builder for model
        solver->init_data(&DS);  // associated with state
        solver->init_calc();
        solver->exec_kernel();
    }
    auto end          = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    // report time cost
    std::cout << solver->scheme(total_time);
    std::cout << "Using total time " << total_time << " s\n";
    return 0;
}

int Handler::run_multiple(Param* P) {
    DataSet DS;

    auto begin = std::chrono::steady_clock::now();
    {
        solver->read_param(P);   // parameter parser and secondary builder for model
        solver->init_data(&DS);  // associated with state
        solver->init_calc(0);

        auto& corr     = Kernel_Record::get_correlation();
        int nframe     = corr.frame;
        int nsize      = corr.size;
        int total_size = corr.size * corr.frame;

        Result corr_sum(corr);
        Result corr_mpi(corr);

        int M = P->get<int>("M", LOC(), 1);
        int istart, iend;
        std::string directory = P->get<std::string>("directory", LOC(), "default");

        MPI_Guard guard{};
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Guard::range(0, M, istart, iend);

        if (MPI_Guard::rank == 0) std::cout << P->repr() << std::endl;

        for (int icycle = istart, icalc = 0; icycle < iend; ++icycle, ++icalc) {
            auto mid1 = std::chrono::steady_clock::now();

            solver->init_calc(icycle);
            solver->exec_kernel(icycle);

            // clear correlation information
            for (int iframe = 0, i = 0; iframe < nframe; ++iframe) {
                bool valid = (corr.stat[iframe] == 1);
                corr_sum.stat[iframe] += valid ? 1 : 0;
                for (int isize = 0; isize < nsize; ++isize, ++i) {
                    corr_sum.data[i] += valid ? corr.data[i] : 0.0e0;
                    corr.data[i] = 0.0e0;
                }
                corr.stat[iframe] = 0;
            }

            auto mid2 = std::chrono::steady_clock::now();
            if (icycle == istart && MPI_Guard::rank == 0) {
                std::cout << "expect time: "
                          << static_cast<std::chrono::duration<double>>(mid2 - mid1).count() * (iend - istart)  //
                          << std::endl;
            }
        }
        MPI_Reduce(corr_sum.stat.data(), corr_mpi.stat.data(), corr.frame, MPI::INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(corr_sum.data.data(), corr_mpi.data.data(), total_size, MPI::DOUBLE_PRECISION, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        for (int iframe = 0, i = 0; iframe < nframe; ++iframe) {
            for (int isize = 0; isize < nsize; ++isize, ++i) corr_sum.data[i] /= corr_sum.stat[iframe];
        }
        corr_sum.save(utils::concat(directory, "/corr-mpi", MPI_Guard::rank, ".dat"), 0, -1, true);
        if (MPI_Guard::isroot) {
            for (int iframe = 0, i = 0; iframe < nframe; ++iframe) {
                for (int isize = 0; isize < nsize; ++isize, ++i) corr_mpi.data[i] /= corr_mpi.stat[iframe];
            }
            corr_mpi.save(utils::concat(directory, "/corr.dat"), 0, -1, true);
        }
    }
    auto end          = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    // report time cost
    if (MPI_Guard::isroot) {
        std::cout << solver->scheme(total_time);
        std::cout << "Using total time " << total_time << " s\n";
    }
    return 0;
}


};  // namespace PROJECT_NS

#include "Handler.h"

#include "core/DataSet.h"
#include "kernels/Kernel_DataSetHandles.h"
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

        auto& corr     = Kernel_Record::get_correlation();
        int total_size = corr.size * corr.frame;

        Result corr_sum(corr);
        Result corr_mpi(corr);

        int N_mc       = P->get<int>("N_mc", LOC(), 1);
        double dt      = P->get<double>("dt", LOC(), phys::time_d);
        double unit_dt = P->get<double>("unit_dt", LOC(), phys::time_d, 1);
        int sstep      = P->get<int>("sstep", LOC(), 1);
        int istart, iend;

        MPI_Guard gaurd{};
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Guard::range(0, N_mc, istart, iend);

        if (MPI_Guard::rank == 0) std::cout << P->repr();

        std::cout << N_mc << ", " << istart << ", " << iend << "\n";
        for (int icycle = istart, icalc = 0; icycle < iend; ++icycle, ++icalc) {
            solver->init_calc(icycle);
            solver->exec_kernel(icycle);

            // clear correlation information
            for (int i = 0; i < total_size; ++i) {
                corr_sum.data[i] += corr.data[i];
                corr.data[i] = 0.0f;
            }
        }
        MPI_Reduce(corr_sum.data.data(), corr_mpi.data.data(), total_size, MPI::DOUBLE_PRECISION, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        for (int i = 0; i < total_size; ++i) corr_sum.data[i] /= (iend - istart);
        corr_sum.save(utils::concat("corr-mpi", MPI_Guard::rank, ".dat"), 0, sstep * dt / unit_dt, true);
        if (MPI_Guard::isroot) {
            for (int i = 0; i < total_size; ++i) corr_mpi.data[i] /= N_mc;
            corr_mpi.save("corr.dat", 0, sstep * dt / unit_dt, true);
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

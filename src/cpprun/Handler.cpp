#include "Handler.h"

#include <chrono>

#include "kids/Kernel_Recorder.h"
#include "kids/ModelFactory.h"
#include "kids/Param.h"
#include "kids/RuleSet.h"
#include "kids/SolverFactory.h"
#include "kids/debug_utils.h"

//
#include "mpi_utils.h"

namespace PROJECT_NS {


Handler::Handler(const std::string& solver_kernel_name, const std::string& model_name) {
    model        = defaultModelFactory(model_name);
    auto solver1 = defaultSolverFactory("Sampling", model);
    auto solver2 = defaultSolverFactory(solver_kernel_name, model);
    solvers.push_back(solver1);
    solvers.push_back(solver2);
    solver1_kernel = solver1->getSolverKernel();
    solver2_kernel = solver2->getSolverKernel();
};

int Handler::run(std::shared_ptr<Param>& PM) {
    std::string hdlr_str = PM->get_string({"solver.handler", "handler"}, LOC(), "");
    if (false) {
    } else if (hdlr_str == "help") {
        run_help(PM);
    } else if (hdlr_str == "help_param") {
        run_help_param(PM);
    } else if (hdlr_str == "help_dataset") {
        run_help_dataset(PM);
    } else if (hdlr_str == "sampling") {
        run_sampling(PM);
    } else if (hdlr_str == "parallel") {
        run_parallel(PM);
    } else if (hdlr_str == "single") {
        // run_single(PM);
    } else if (hdlr_str == "single_mpi") {
        // run_single_mpi(PM);
    } else {
        throw std::runtime_error("unknown handler type!");
    }
    return 0;
}

int Handler::run_single(std::shared_ptr<Param>& PM) {
    std::shared_ptr<DataSet> DS    = std::shared_ptr<DataSet>(new DataSet());
    auto                     begin = std::chrono::steady_clock::now();
    Status                   stat;
    {
        solver2_kernel->setInputParam(PM);
        solver2_kernel->setInputDataSet(DS);
        solver2_kernel->initializeKernel(stat);
        solver2_kernel->executeKernel(stat);
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    std::cout << solver2_kernel->generateInformationString(total_time);
    std::cout << "Using total time " << total_time << " s\n";
    return 0;
}

int Handler::run_single_mpi(std::shared_ptr<Param>& PM) {
    std::shared_ptr<DataSet> DS = std::shared_ptr<DataSet>(new DataSet());
    Status                   stat;
    auto                     begin = std::chrono::steady_clock::now();
    {
        MPI_Guard guard(1);
        MPI_Barrier(MPI_COMM_WORLD);

        solver2_kernel->setInputParam(PM);
        solver2_kernel->setInputDataSet(DS);
        solver2_kernel->initializeKernel(stat);
        solver2_kernel->executeKernel(stat);
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    std::cout << solver2_kernel->generateInformationString(total_time);
    std::cout << "Using total time " << total_time << " s\n";
    return 0;
}

int Handler::run_parallel(std::shared_ptr<Param>& PM) {
    std::shared_ptr<DataSet> DS = std::shared_ptr<DataSet>(new DataSet());
    // std::vector<shared_ptr<RuleSet>> RSL;
    // RSL.push(std::shared_ptr<RuleSet>(new RuleSet()));

    // bool  not_parsed = _ruleset->getRules().size() == 0;
    // auto& json       = *(_param->dataObject());
    // if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
    //     for (auto& j : (json["result"])) token(j);
    // }

    std::cout << solver1_kernel->generateInformationString(1.0);
    std::cout << solver2_kernel->generateInformationString(1.0);

    auto   begin = std::chrono::steady_clock::now();
    Status stat;
    {
        solver1_kernel->setInputParam(PM);
        solver2_kernel->setInputParam(PM);
        solver1_kernel->setInputDataSet(DS);
        solver2_kernel->setInputDataSet(DS);

        // get Monte Carlo Dimension
        MPI_Guard guard(solver1_kernel->montecarlo);
        MPI_Barrier(MPI_COMM_WORLD);

        if (MPI_Guard::isroot) std::cout << PM->repr() << std::endl;

        for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
            auto mid1 = std::chrono::steady_clock::now();

            stat.icalc = icalc;

            solver1_kernel->initializeKernel(stat);
            solver1_kernel->executeKernel(stat);
            solver1_kernel->finalizeKernel(stat);

            solver2_kernel->initializeKernel(stat);
            solver2_kernel->executeKernel(stat);
            solver2_kernel->finalizeKernel(stat);

            auto mid2 = std::chrono::steady_clock::now();
            if (icalc == guard.istart && MPI_Guard::isroot) {
                std::cout << "expect time: "
                          << static_cast<std::chrono::duration<double>>(mid2 - mid1).count() *
                                 (guard.iend - guard.istart)  //
                          << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // @todo customizable collection
        auto collect = solver2_kernel->getRuleSet()->getCollect().data();
        auto reduced = solver2_kernel->getRuleSet()->getReduced().data();
        // @bad because it should in public domain, but collect return null for blank mpi
        for (int i = 0; i < collect.size(); ++i) {
            // std::cout << std::get<0>(collect[i]) << "\n";
            // std::cout << std::get<0>(reduced[i]) << "\n";
            auto [key1, from_data, type1, size1, nframe1] = collect[i];
            auto [key2, to_data, type2, size2, nframe2]   = reduced[i];
            MPI_Guard::reduce(std::make_tuple(type1, from_data, to_data, size1));
        }
        // report time cost
        if (MPI_Guard::isroot) { RuleSet::flush_all(solver2_kernel->directory, 2); }
        // std::cout << DS->repr() << "\n";
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    // report time cost
    if (MPI_Guard::isroot) {
        std::cout << solver2_kernel->generateInformationString(total_time);
        std::cout << "Using total time " << total_time << " s\n";
    }
    return 0;
}

int Handler::run_sampling(std::shared_ptr<Param>& PM) {
    std::shared_ptr<DataSet> DS    = std::shared_ptr<DataSet>(new DataSet());
    auto                     begin = std::chrono::steady_clock::now();
    Status                   stat;
    {
        solver1_kernel->setInputParam(PM);
        solver1_kernel->setInputDataSet(DS);
        solver1_kernel->initializeKernel(stat);

        MPI_Guard guard(solver1_kernel->montecarlo);
        MPI_Barrier(MPI_COMM_WORLD);

        if (MPI_Guard::rank == 0) std::cout << PM->repr() << std::endl;

        for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
            auto mid1  = std::chrono::steady_clock::now();
            stat.icalc = icalc;
            solver1_kernel->initializeKernel(stat);
        }
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    // report time cost
    if (MPI_Guard::isroot) {
        std::cout << solver1_kernel->generateInformationString(total_time);
        std::cout << "Using total time " << total_time << " s\n";
    }
    return 0;
}

int Handler::run_help_dataset(std::shared_ptr<Param>& PM) {
    {
        for (auto&& i : VARIABLE_BASE::_LIST) {  //
            std::cout << i->name() << ", " << i->doc() << "\n";
        }
    }
    return 0;
}


};  // namespace PROJECT_NS

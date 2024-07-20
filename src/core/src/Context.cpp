#include "kids/Context.h"

#include <chrono>

namespace PROJECT_NS {

Context::Context(std::shared_ptr<Platform> plat, std::shared_ptr<System> sys,
                 std::vector<std::shared_ptr<Solver>> solvers)
    : _platform{plat}, _system{sys}, _solvers{solvers} {
    //
}

Status& Context::run_all(Status& stat) {
    while (stat.istage < _solvers.size()) run(stat);
    return stat;
}

Status& Context::run(Status& stat) {
    auto& solvers = _solvers[stat.istage];
    auto& PM      = _param;
    auto& DS      = _dataset;
    auto  begin   = std::chrono::steady_clock::now();
    // {
    //     // initialization of Param & DataSet for solvers
    //     for (auto& solver : solvers) {
    //         auto solver_kernel = solver->getSolverKernel();
    //         solver_kernel->setInputParam(PM);
    //         solver_kernel->setInputDataSet(DS);
    //     }

    //     // get Monte Carlo Dimension by platform
    //     // auto&& guard = _platform.generateGaurd();
    //     // guard.barrier();
    //     // if (guard.isroot) std::cout << PM->repr() << std::endl;

    //     MPI_Guard guard(solver_kernel->montecarlo);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (MPI_Guard::isroot) std::cout << PM->repr() << std::endl;

    //     // execute for solver_kernels
    //     for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
    //         stat.icalc = icalc;

    //         int isolver = 0;
    //         for (auto& solver : solvers) {
    //             stat.isolver = isolver;

    //             auto& solver_kernel = solver->getSolverKernel();
    //             solver_kernel->initializeKernel(stat);
    //             solver_kernel->executeKernel(stat);
    //             solver_kernel->finalizeKernel(stat);
    //             ++isolver;
    //         }
    //     }

    //     MPI_Barrier(MPI_COMM_WORLD);
    //     std::cout << guard.istart << ";" << guard.iend << " !\n";

    //     // initialization of Param & DataSet for solvers
    //     for (auto& solver : solvers) {
    //         auto solver_kernel = solver->getSolverKernel();
    //         auto collect       = solver_kernel->getRuleSet()->getCollect().data();
    //         auto reduced       = solver_kernel->getRuleSet()->getReduced().data();
    //         // @bad because it should in public domain, but collect return null for blank mpi
    //         for (int i = 0; i < collect.size(); ++i) {
    //             std::cout << std::get<0>(collect[i]) << "\n";
    //             std::cout << std::get<0>(reduced[i]) << "\n";
    //             auto [key1, from_data, type1, size1, nframe1] = collect[i];
    //             auto [key2, to_data, type2, size2, nframe2]   = reduced[i];
    //             MPI_Guard::reduce(std::make_tuple(type1, from_data, to_data, size1));
    //         }
    //     }

    //     // report time cost
    //     if (MPI_Guard::isroot) { RuleSet::flush_all(solver_kernel->directory, 2); }
    // }
    // auto   end        = std::chrono::steady_clock::now();
    // double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();

    // // report time cost
    // if (MPI_Guard::isroot) {
    //     for (auto& solver : solvers) {
    //         auto solver_kernel = solver->getSolverKernel();
    //         std::cout << solver_kernel->generateInformationString(total_time);
    //         std::cout << "Using total time " << total_time << " s\n";
    //     }
    // }

    stat.istage++;
    return stat;
}


};  // namespace PROJECT_NS

#include "kids/Context.h"

namespace PROJECT_NS {

Context::Context(std::shared_ptr<Platform> plat, std::shared_ptr<System> sys,
                 std::vector<std::shared_ptr<Solver>> solvers)
    : _platform{plat}, _system{sys}, _solvers{solvers} {
    //
}

Status& Context::execute(Status& stat) {
    // auto   begin = std::chrono::steady_clock::now();
    // {
    //     for (auto& solver : _solvers) solver->setInputParam(_param);
    //     for (auto& solver : _solvers) solver->setInputDataSet(_dataset);
    //     solver->initializeKernel(stat);  // @necessary?

    //     // get Monte Carlo Dimension from Param
    //     int         M         = PM->get_int("M", LOC(), 1);
    //     std::string directory = PM->get_string("directory", LOC(), "default");
    //     MPI_Guard   guard(M);
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
    //         stat.icalc = icalc;
    //         solver->initializeKernel(stat);
    //         solver->executeKernel(stat);
    //         solver->finalizeKernel(stat);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);

    //     auto collect = solver->getRuleSet()->getResult(1).data();
    //     auto reduced = solver->getRuleSet()->getResult(2).data();
    //     for (int i = 0; i < collect.size(); ++i) {
    //         std::cout << std::get<0>(collect[i]) << "\n";
    //         std::cout << std::get<0>(reduced[i]) << "\n";
    //         auto [key1, from_data, type1, size1, nframe1] = collect[i];
    //         auto [key2, to_data, type2, size2, nframe2]   = reduced[i];
    //         MPI_Guard::reduce(std::make_tuple(type1, from_data, to_data, size1));
    //     }
    //     // report time cost
    //     if (MPI_Guard::isroot) RuleSet::flush_all(directory, 2);
    // }
    // auto   end        = std::chrono::steady_clock::now();
    // double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();
    return stat;
}


};  // namespace PROJECT_NS

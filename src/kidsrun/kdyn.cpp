#include <gflags/gflags.h>
#include <glog/logging.h>
// #include <gperftools/profiler.h>

#include <iomanip>
#include <iostream>

#include "Handler.h"
#include "defined_gflags.h"
#include "ghc/filesystem.hpp"
#include "kids/Kernel.h"
#include "kids/Model.h"
#include "kids/Param.h"
#include "kids/Solver.h"
#include "version.h"

using namespace PROJECT_NS;
namespace fs = ghc::filesystem;

int main(int argc, char* argv[]) {
    /* profiling settings */
    // plSetFilename("record.pltraw");
    // plInitAndStart("opendf", PL_MODE_STORE_IN_FILE);
    // plDeclareThread("Main");

    /* setup gflags & glog */
    gflags::SetVersionString(repo_version);
    gflags::SetUsageMessage("Kernel Integrated Dynamics Simulator");
    gflags::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_INFO);
    google::SetLogDestination(google::GLOG_INFO, utils::concat("./", FLAGS_d, "/").c_str());
    google::SetLogFilenameExtension(".log");
    // FLAGS_timestamp_in_logfile_name = false;
    // FLAGS_log_prefix = false;
    // FLAGS_colorlogtostderr          = true;
    FLAGS_logtostderr               = 0;
    FLAGS_alsologtostderr           = 0;
    FLAGS_logbufsecs                = 0;
    FLAGS_max_log_size              = 5;
    FLAGS_stop_logging_if_full_disk = true;

    /* KIDS simulation */
    Status stat;

    std::shared_ptr<Param> PM = std::shared_ptr<Param>(new Param(FLAGS_p, Param::fromFile));
    check_and_sync_from_flags(PM);
    std::string model_name  = PM->get_string({"model.name"}, LOC());
    std::string solver_name = PM->get_string({"solver.name"}, LOC());

    std::shared_ptr<Model>   model         = defaultModelFactory(model_name);
    std::shared_ptr<Solver>  solver        = defaultSolverFactory(solver_name, model);
    std::shared_ptr<Kernel>  solver_kernel = solver->getSolverKernel();
    std::shared_ptr<DataSet> DS            = std::shared_ptr<DataSet>(new DataSet());

    auto begin = std::chrono::steady_clock::now();
    {
        solver_kernel->setInputParam(PM);
        solver_kernel->setInputDataSet(DS);
        solver_kernel->initializeKernel(stat);
        Simple_Guard guard(solver_kernel->montecarlo);
        for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
            solver_kernel->initializeKernel(stat);
            solver_kernel->executeKernel(stat);
            solver_kernel->finalizeKernel(stat);
        }
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();
    std::cout << solver_kernel->generateInformationString(total_time);
    std::cout << "Using total time " << total_time << " s\n";

    /* finalize */
    gflags::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();
    // ProfilerStop();
    // plStopAndUninit();

    return 0;
}

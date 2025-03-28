#include <iomanip>
#include <iostream>

#include "defined_gflags.h"
#include "kids/Kernel.h"
#include "kids/Model.h"
#include "kids/ModelFactory.h"
#include "kids/Param.h"
#include "kids/Solver.h"
#include "kids/SolverFactory.h"
#include "simple_guard.h"
#include "version.h"

using namespace PROJECT_NS;

int main(int argc, char* argv[]) {
    ////////////////////////////////////////////////////////////////////////////////////////
    /* profiling settings */
    // plSetFilename("record.pltraw");
    // plInitAndStart("opendf", PL_MODE_STORE_IN_FILE);
    // plDeclareThread("Main");
    /* gflags & glog settings */
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

    ////////////////////////////////////////////////////////////////////////////////////////

    /* KIDS simulation */
    Status stat;

    std::shared_ptr<Param> PM = std::shared_ptr<Param>(new Param(FLAGS_p, Param::fromFile));
    check_and_sync_from_gflags(PM);
    PM->set_string("dump", "samp");
    PM->set_string("solver.dump", "samp");
    std::ofstream ofs(utils::concat(FLAGS_d, "/", "input"));
    ofs << PM->repr();
    ofs.close();

    std::string model_name = PM->get_string({"model.name"}, LOC());
    std::size_t BGIDX      = PM->get_int({"BGIDX"}, LOC(), 0);

    std::shared_ptr<Model>   model         = defaultModelFactory(model_name);
    std::shared_ptr<Solver>  solver        = defaultSolverFactory("Sampling", model);
    std::shared_ptr<Kernel>  solver_kernel = solver->getSolverKernel();
    std::shared_ptr<DataSet> DS            = std::shared_ptr<DataSet>(new DataSet());

    auto begin = std::chrono::steady_clock::now();
    {
        solver_kernel->setInputParam(PM);
        solver_kernel->setInputDataSet(DS);
        Simple_Guard guard(BGIDX, solver_kernel->montecarlo);
        for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
            stat.icalc = icalc;
            solver_kernel->initializeKernel(stat);
            solver_kernel->executeKernel(stat);
            solver_kernel->finalizeKernel(stat);
        }
    }
    auto   end        = std::chrono::steady_clock::now();
    double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();
    std::cout << solver_kernel->generateInformationString(total_time);
    std::cout << "Using total time " << total_time << " s\n";

    ////////////////////////////////////////////////////////////////////////////////////////
    /* finalize */
    gflags::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();
    // ProfilerStop();
    // plStopAndUninit();
    return 0;
}

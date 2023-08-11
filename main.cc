#include <gflags/gflags.h>
#include <glog/logging.h>
// #include <gperftools/profiler.h>

#include <iomanip>
#include <iostream>

#include "Handler.h"
#include "core/Param.h"
#include "generate/version.h"
#include "thirdpart/ghc/filesystem.hpp"

DEFINE_string(p, "param.json", "paramemter inputs");
DEFINE_string(d, "default", "directory for simulation");
DEFINE_string(s, "ENER,TRAJ,ETRAJ,CORR", "outputers; seperated by comma");
DEFINE_bool(w, false, "rewrite output");
DEFINE_bool(timing, false, "timing kernels");
DEFINE_bool(trace, false, "trace into outfilestream");

DEFINE_int32(nsave_mpi, 5, "save backups during monte carlo");
DEFINE_int32(nsave_time, 5, "save backups during a simulation");
DEFINE_double(everysave, 1.0f, "interval for saving a context (unit in hours)");

using namespace opendf;
namespace fs = ghc::filesystem;

int main(int argc, char* argv[]) {
    /* profiling settings */
    // ProfilerStart("demo.prof");

    // plSetFilename("record.pltraw");
    // plInitAndStart("opendf", PL_MODE_STORE_IN_FILE);
    // plDeclareThread("Main");

    /* setup gflags */
    gflags::SetVersionString(repo_version);
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    /* read parameter file (json format) */
    Param P(FLAGS_p, Param::fromFile);

    /* creat directory for simulation */
    FLAGS_d = P.get<std::string>("jobid", LOC(), FLAGS_d);
    if (fs::exists(FLAGS_d) && FLAGS_w == false) {
        throw std::runtime_error(
            utils::concat("Working directory = [", FLAGS_d, "] already exists. Please specify -w to force start.\n"));
    }
    fs::create_directory(FLAGS_d);

    auto&& j = *(P.pjson());
    if (j.count("model_file") > 0 && j.count("model_id") > 0) {
        Param TEMP(j["model_file"].as_string(), Param::fromFile);
        j["model_param"] = (*(TEMP.pjson()))[j["model_id"].as_string()];
    }
    if (j.count("solver_file") > 0 && j.count("solver_id") > 0) {
        Param TEMP(j["solver_file"].as_string(), Param::fromFile);
        j["solver_param"] = (*(TEMP.pjson()))[j["solver_id"].as_string()];
    }

    /* setup glog */
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::GLOG_INFO);
    google::SetLogDestination(google::GLOG_INFO, utils::concat("./", FLAGS_d, "/").c_str());
    google::SetLogFilenameExtension(".log");
    // FLAGS_timestamp_in_logfile_name = false;
    // FLAGS_log_prefix = false;
    // FLAGS_colorlogtostderr          = true;  // Set log color
    FLAGS_logtostderr               = 0;
    FLAGS_alsologtostderr           = 0;
    FLAGS_logbufsecs                = 0;     // Set log output speed(s)
    FLAGS_max_log_size              = 5;     // Set max log file size
    FLAGS_stop_logging_if_full_disk = true;  // If disk is full

    /* task block */
    std::string model_name   = P.get<std::string>("model", LOC());
    std::string solver_name  = P.get<std::string>("solver", LOC());
    std::string handler_name = P.get<std::string>("handler", LOC(), "single");

    (*P.pjson())["trace"]     = FLAGS_trace;
    (*P.pjson())["is_timing"] = FLAGS_timing;

    Handler myhandler = Handler(Handler::multiple, solver_name, model_name);

    myhandler.run(&P);

    /* finalize */
    gflags::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();

    // ProfilerStop();
    // plStopAndUninit();

    return 0;
}

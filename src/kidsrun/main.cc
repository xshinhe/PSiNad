#include <gflags/gflags.h>
#include <glog/logging.h>
// #include <gperftools/profiler.h>

#include <iomanip>
#include <iostream>

#include "Handler.h"
#include "ghc/filesystem.hpp"
#include "kids/Param.h"
#include "version.h"

DEFINE_bool(w, false, "Enables rewritting the output");
DEFINE_string(p, "param.json", "paramemter inputs");

DEFINE_string(handler,     //
              "parallel",  //
              "Specifies the handler type"
              "\n[parallel | single | single_mpi | sampling | help | help_param | help_dataset ]");
DEFINE_string(d, "default", "Specifies the output directory");
DEFINE_string(load, "", "Specifies the dataset file to load");
DEFINE_string(dump, "", "Specifies the dataset file to dump");
DEFINE_double(backup_time, -1.0, "Specifies the timestep for backup (/1h)");
DEFINE_bool(timing, false, "Enables profiling for time costs");

using namespace PROJECT_NS;
namespace fs = ghc::filesystem;

int main(int argc, char* argv[]) {
    /* profiling settings */
    // ProfilerStart("demo.prof");

    // plSetFilename("record.pltraw");
    // plInitAndStart("opendf", PL_MODE_STORE_IN_FILE);
    // plDeclareThread("Main");

    /* setup gflags */
    gflags::SetVersionString(repo_version);
    gflags::SetUsageMessage("Kernel Integrated Dynamics Simulator");
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    /* read parameter file (json format) */
    Param  PM(FLAGS_p, Param::fromFile);
    auto&& j = *(PM.pjson());

    j["directory"]   = FLAGS_d;
    j["timing"]      = FLAGS_timing;
    j["handler"]     = FLAGS_handler;
    j["backup_time"] = FLAGS_backup_time;
    if (FLAGS_load != "") j["load"] = FLAGS_load;
    if (FLAGS_dump != "") j["dump"] = FLAGS_dump;

    /* creat directory for simulation */
    if (fs::exists(FLAGS_d) && FLAGS_w == false) {
        throw std::runtime_error(
            utils::concat("Working directory = [", FLAGS_d, "] already exists. Please specify -w to force start.\n"));
    }
    try {
        fs::create_directory(FLAGS_d);  // sometime it raise bugs
    } catch (std::runtime_error& e) {
        throw std::runtime_error("create_directory failed");
        std::cout << "some error!!!\n";
    }

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
    std::string model_name   = PM.get_string("model", LOC());
    std::string solver_name  = PM.get_string("solver", LOC());
    std::string handler_name = PM.get_string("handler", LOC(), "single");

    Handler myhandler = Handler(solver_name, model_name);
    myhandler.run(&PM);

    /* finalize */
    gflags::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();

    // ProfilerStop();
    // plStopAndUninit();

    return 0;
}

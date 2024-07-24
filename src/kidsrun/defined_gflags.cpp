#include "defined_gflags.h"

#include "ghc/filesystem.hpp"
#include "kids/concat.h"

using namespace PROJECT_NS;
namespace fs = ghc::filesystem;

DEFINE_bool(w, false, "Enables rewritting the output");
DEFINE_string(p, "param.json", "paramemter inputs");

DEFINE_string(handler,     //
              "parallel",  //
              "Specifies the handler type\n"
              "[parallel | single | single_mpi | sampling | help | help_param | help_dataset ]");
DEFINE_string(d, "default", "Specifies the output directory");
DEFINE_string(load, "", "Specifies the dataset file to load");
DEFINE_string(dump, "", "Specifies the dataset file to dump");
DEFINE_double(backup_time, -1.0, "Specifies the timestep for backup (/1h)");
DEFINE_bool(timing, false, "Enables simple profiling for time costs");
DEFINE_bool(profiling, false, "Enables high-performance profiling for time costs");
DEFINE_int32(seed, -1, "random seed");
DEFINE_int32(BGIDX, 0, "traj idx start");

void check_and_sync_from_gflags(std::shared_ptr<Param> PM) {
    PM->set_string("directory", FLAGS_d);
    PM->set_bool("timing", FLAGS_timing);
    PM->set_string("handler", FLAGS_handler);
    PM->set_real("backup_time", FLAGS_backup_time);
    if (FLAGS_load != "") PM->set_string("load", FLAGS_load);
    if (FLAGS_dump != "") PM->set_string("dump", FLAGS_dump);
    if (FLAGS_seed != -1) PM->set_int("seed", FLAGS_seed);
    if (FLAGS_BGIDX != 0) PM->set_int("BGIDX", FLAGS_BGIDX);

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
}

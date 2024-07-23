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


void check_and_sync_from_gflags(std::shared_ptr<Param> PM) {
    auto&& j         = *(PM->pjson());
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
}

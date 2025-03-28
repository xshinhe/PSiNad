#include <iomanip>
#include <iostream>

#include "defined_gflags.h"
#include "kids/Kernel.h"
#include "kids/Model.h"
#include "kids/ModelFactory.h"
#include "kids/Param.h"
#include "kids/Solver.h"
#include "kids/SolverFactory.h"
#include "mpi_guard.h"
#include "version.h"

using namespace PROJECT_NS;

int main(int argc, char* argv[])
{
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
    google::SetLogDestination(google::GLOG_INFO,
                              utils::concat("./", FLAGS_d, "/").c_str());
    google::SetLogFilenameExtension(".log");
    // FLAGS_timestamp_in_logfile_name = false;
    // FLAGS_log_prefix = false;
    // FLAGS_colorlogtostderr          = true;
    FLAGS_logtostderr = 0;
    FLAGS_alsologtostderr = 0;
    FLAGS_logbufsecs = 0;
    FLAGS_max_log_size = 5;
    FLAGS_stop_logging_if_full_disk = true;

    ////////////////////////////////////////////////////////////////////////////////////////

    /* KIDS simulation */
    Status stat;

    std::shared_ptr<Param> PM =
        std::shared_ptr<Param>(new Param(FLAGS_p, Param::fromFile));
    check_and_sync_from_gflags(PM);

    std::string model_name = PM->get_string({"model.name"}, LOC());
    std::string solver_name = PM->get_string({"solver.name"}, LOC());
    std::string solver_scheme = PM->get_string({"solver.scheme"}, LOC(), "");
    std::size_t BGIDX = PM->get_int({"BGIDX"}, LOC(), 0);

    // apply_scheme(PM, solver_scheme);
    std::ofstream ofs(utils::concat(FLAGS_d, "/", "input"));
    ofs << PM->repr();
    ofs.close();

    std::shared_ptr<Model> model = defaultModelFactory(model_name);
    std::shared_ptr<Solver> solver1 = defaultSolverFactory("Sampling", model);
    std::shared_ptr<Kernel> solver1_kernel = solver1->getSolverKernel();
    std::shared_ptr<Solver> solver2 = defaultSolverFactory(solver_name, model);
    std::shared_ptr<Kernel> solver2_kernel = solver2->getSolverKernel();
    std::shared_ptr<DataSet> DS = std::shared_ptr<DataSet>(new DataSet());

    auto begin = std::chrono::steady_clock::now();
    {
        solver1_kernel->setInputParam(PM);
        solver2_kernel->setInputParam(PM);
        solver1_kernel->setInputDataSet(DS);
        solver2_kernel->setInputDataSet(DS);

        MPI_Guard guard(BGIDX, solver1_kernel->montecarlo);

        if (MPI_Guard::isroot)
        {
            std::cout << solver1_kernel->generateInformationString(1.0);
            std::cout << solver2_kernel->generateInformationString(1.0);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        for (int icalc = guard.istart; icalc < guard.iend; ++icalc)
        {
            auto mid1 = std::chrono::steady_clock::now();

            stat.icalc = icalc;

            solver1_kernel->initializeKernel(stat);
            solver1_kernel->executeKernel(stat);
            solver1_kernel->finalizeKernel(stat);

            solver2_kernel->initializeKernel(stat);
            solver2_kernel->executeKernel(stat);
            solver2_kernel->finalizeKernel(stat);

            int flevel = -1;
            if (stat.icalc == guard.istart) flevel = -2;
            RuleSet::flush_all(solver2_kernel->directory,
                               utils::concat(".mpi", MPI_Guard::rank), flevel);

            auto mid2 = std::chrono::steady_clock::now();
            if (icalc == guard.istart && MPI_Guard::isroot)
            {
                std::cout << "expect time: "
                          << static_cast<std::chrono::duration<double>>(mid2 -
                                                                        mid1)
                                     .count() *
                                 (guard.iend - guard.istart) //
                          << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        auto collect = solver2_kernel->getRuleSet()->getCollect().data();
        auto reduced = solver2_kernel->getRuleSet()->getReduced().data();
        for (int i = 0; i < collect.size(); ++i)
        {
            auto [key1, from_data, type1, size1, nframe1] = collect[i];
            auto [key2, to_data, type2, size2, nframe2] = reduced[i];
            MPI_Guard::reduce(
                std::make_tuple(type1, from_data, to_data, size1));
        }
        // RuleSet::flush_all(solver2_kernel->directory, utils::concat(".mpi",
        // MPI_Guard::rank), 1);
        if (MPI_Guard::isroot)
        {
            RuleSet::flush_all(solver2_kernel->directory, "", 2);
        }
    }
    auto end = std::chrono::steady_clock::now();
    double total_time =
        static_cast<std::chrono::duration<double>>(end - begin).count();

    if (MPI_Guard::isroot)
    {
        std::cout << solver2_kernel->generateInformationString(total_time);
        std::cout << "Using total time " << total_time << " s\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    /* finalize */
    gflags::ShutDownCommandLineFlags();
    google::ShutdownGoogleLogging();
    // ProfilerStop();
    // plStopAndUninit();
    return 0;
}

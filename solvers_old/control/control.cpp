#include "control.h"

#include "../utils/definitions.h"
#include "../utils/hdf5_utils.h"

Control::Control(Solver* pS) : pSolver(pS) {
    // refer random seed state
    if (FLAGS_r == "") {
        pContext = new File("run.h5", File::ReadWrite | File::Create | File::Truncate);
        rand_rng = SeededEngine(nullptr);
    } else {
        if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "run.h5");
        try {
            pContext = new File("run.h5", File::ReadWrite);
        } catch (HighFive::Exception& e) { LOG(FATAL); }
        if (FLAGS_read_seed) {
            LOAD(pContext, rng_seeds, rng_t::state_size);
            rand_rng = SeededEngine(rng_seeds);
        } else {
            rand_rng = SeededEngine(nullptr);
        }
    }
};

Control::~Control(){
    // save random seed state
};

SingleControl::SingleControl(Solver* pS) : Control(pS){};

int SingleControl::run() {
    if (FLAGS_r == "") {
        pContext = new File("run.h5", File::ReadWrite | File::Create | File::Truncate);
        rand_rng = SeededEngine(nullptr);
    } else {
        if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "run.h5");
        try {
            pContext = new File("run.h5", File::ReadWrite);
        } catch (HighFive::Exception& e) { LOG(FATAL); }
        if (FLAGS_read_seed) {
            LOAD(pContext, rng_seeds, rng_t::state_size);
            rand_rng = SeededEngine(rng_seeds);
        } else {
            rand_rng = SeededEngine(nullptr);
        }
    }

    pSolver->run_impl();

    // save random seeds
    DUMP(pContext, rng_seeds, rng_t::state_size);

    pContext->flush();
    delete pContext;  // bugs?

    return 0;
}


MPIControl::MPIControl(Solver* pS) : Control(pS){};

int MPIControl::run() {
    mpi_utils_init();
    if (FLAGS_r == "") {
        pContext = new File("run.h5", File::ReadWrite | File::Create | File::Truncate,
                            MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

        rand_rng = SeededEngine(nullptr);
    } else {
        if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "run.h5");
        try {
            pContext = new File("run.h5", File::ReadWrite, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
        } catch (HighFive::Exception& e) { LOG(FATAL); }

        if (FLAGS_read_seed) {
            LOAD(pContext, rng_seeds, rng_t::state_size);
            rand_rng = SeededEngine(rng_seeds);
        } else {
            rand_rng = SeededEngine(nullptr);
        }
    }

    pSolver->run_impl();

    // save random seeds
    DUMP(pContext, rng_seeds, rng_t::state_size);

    pContext->flush();
    delete pContext;  // bugs?
    mpi_utils_final();

    return 0;
}


MonteCarloControl::MonteCarloControl(Solver* pS) : Control(pS){};
int MonteCarloControl::run() {
    mpi_utils_init();
    if (FLAGS_r == "") {
        pContext = new File("run.h5", File::ReadWrite | File::Create | File::Truncate,
                            MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

        rand_rng = SeededEngine(nullptr);
    } else {
        if (mpi_rank == 0) utils::copyfile_from_to(FLAGS_r, "run.h5");
        try {
            pContext = new File("run.h5", File::ReadWrite, MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));
        } catch (HighFive::Exception& e) { LOG(FATAL); }

        if (FLAGS_read_seed) {
            LOAD(pContext, rng_seeds, rng_t::state_size);
            rand_rng = SeededEngine(rng_seeds);
        } else {
            rand_rng = SeededEngine(nullptr);
        }
    }

    int nsave = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    if (nsave > FLAGS_nsave_mpi) nsave = FLAGS_nsave_mpi;
    num_real sampunit = dt * sstep * iou.time;

    pSolver->init(-1);  // get occ0

    NAD_TCFer coll    = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collsum = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collmpi = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);

    for (int isave = 0; isave < nsave; ++isave) {
        int eachstart = (isave * ntraj) / nsave, eachend = ((isave + 1) * ntraj) / nsave, istart, iend;

        mpi_range(eachstart, eachend, mpi_nprocs, mpi_rank, istart, iend);
        MPI_Barrier(MPI_COMM_WORLD);
        for (int icycle = istart; icycle < iend; ++icycle) {
            pSolver->itraj = icycle;
            coll.Clear();
            pSolver->run_impl();
            traj(coll, N, F);
            collsum.Amount(coll);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        collmpi.MPIAmount(collsum);

        // MPI_Reduce(Collsum, Collmpi, ncoll, MPI::DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI_Reduce(Statsum, Statmpi, nsamp, MPI::INT, MPI_SUM, 0, MPI_COMM_WORLD);

        if (mpi_rank == 0) collmpi.report(save + "-cache" + std::to_string(isave), sampunit);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    collmpi.MPIAmount(collsum);
    if (mpi_rank == 0 && global::ofs_is_open(global::OFS::CORR)) collmpi.report(save, sampunit);
    collsum.report(save + "-mpi" + std::to_string(mpi_rank), sampunit);
    final(-1);
    return 0;

    // save random seeds
    DUMP(pContext, rng_seeds, rng_t::state_size);

    pContext->flush();
    delete pContext;
    mpi_utils_final();
}

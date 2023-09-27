#include "traj.h"

#include "../utils/definitions.h"

#define if_debug false  // to locate bug
using TCFnucl = num_real;

Traj_Solver::Traj_Solver(const Param& iparm, Model* pM) : Solver(iparm, pM) {
    // cast parameters

    pForceField = dynamic_cast<BO_ForceField*>(pM);

    Param_GetV(level, iparm, 1);    // default level (force level)
    Param_GetV(ntraj, iparm, 1);    // no. of trajectory
    Param_GetV(nstep, iparm, 1);    // no. of steps of each trajectory
    Param_GetV(sstep, iparm, 1);    // skip steps for each sampling
    Param_GetV(nsamp, iparm, 1);    // skip steps for each sampling
    Param_GetV(nrespa, iparm, 10);  // for multiple time-scaling
    Param_GetV(dt, iparm, -1.0f);   // normal time-step
    Param_GetV(tend, iparm, 1.0f);  // total time for simulation

    // update (dt, tend) parameters
    if (pForceField->Suggest_dt() > 0.0f) {
        dt = pForceField->Suggest_dt();  // suggested value from forcefield
    }
    if (pForceField->Suggest_tend() > 0.0f) {
        tend = pForceField->Suggest_tend();  // suggested value from forcefield
    }
    Param_Reset(nstep, sstep * (int(tend / (sstep * dt)) + 1));
    Param_Reset(nsamp, nstep / sstep + 1);
    Param_Reset(halfdt, 0.5f * dt);
    Param_Reset(smalldt, dt / nrespa);
    Param_Reset(largedt, nrespa * dt);

    // get dimensions
    N     = pForceField->get_N();  // get no. of nuclear DOFs
    NN    = N * N;                 // NN := N * N
    nspec = pForceField->nspec();  // get no. of specified dividings
    Ndim  = pForceField->get_Ndim();
    Natom = N / Ndim;

    // Thermostat
    pThermo = new Thermostat(iparm);
    pThermo->init_alloc(N);  // for basic MD, only need N dofs
    try {
        /* @semi-classical: not necessary. And it can output to ofs_ENER
        ALLOCATE_PTR_TO_VECTOR(Htot, nsamp);
        ALLOCATE_PTR_TO_VECTOR(Ltot, nsamp);
        ALLOCATE_PTR_TO_VECTOR(Stot, nsamp);
        ALLOCATE_PTR_TO_VECTOR(Ktot, nsamp);
        ALLOCATE_PTR_TO_VECTOR(Vtot, nsamp);
        */
        ///< sampling variables
        ALLOCATE_PTR_TO_VECTOR(nr0, N);
        ALLOCATE_PTR_TO_VECTOR(np0, N);
        ///< dynamical variables
        ALLOCATE_PTR_TO_VECTOR(nr, N);
        ALLOCATE_PTR_TO_VECTOR(np, N);
        ALLOCATE_PTR_TO_VECTOR(nm, N);
        ALLOCATE_PTR_TO_VECTOR(nf, N);
        ///< forcefield variables
        ALLOCATE_PTR_TO_VECTOR(vpes, 1);
        ALLOCATE_PTR_TO_VECTOR(grad, N);
        ALLOCATE_PTR_TO_VECTOR(hess, NN);

        ALLOCATE_PTR_TO_VECTOR(workr, NN);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    CHECK_GT(dt, 0);
    CHECK_GT(tend, 0);
    CHECK_GT(sstep, 0);
    CHECK_GT(ntraj, 0);
    CHECK_GT(nstep, 0);
    CHECK_GT(nsamp, 0);
    LOG_IF(INFO, if_debug) << "solver initialized";
}

Traj_Solver::~Traj_Solver() { delete pThermo; };

int Traj_Solver::ff_calc1(const int& level) {
    int succ = pForceField->ForceField_npes(vpes, grad, hess, nr, np, level, N, itraj, istep);
    LOG_IF(INFO, if_debug) << "forcefield calculation done";
    for (int i = 0; i < N; ++i) nf[i] = grad[i];
    return succ;
}

int Traj_Solver::init_ofs(const int& itraj) {
    if (itraj < 0) {
        if (global::ofs_is_open(global::OFS::SAMP)) {
            auto fname = utils::concat(save, "-", mpi_rank, ".samp");
            utils::clearFile(fname);
            ofs_SAMP.open(fname);
        }
    } else {  // initialization for each trajectory
        if (global::ofs_is_open(global::OFS::ENER)) {
            auto fname = utils::concat(save, "-", itraj, ".ener");
            utils::clearFile(fname);
            ofs_ENER.open(fname);
        }
        if (global::ofs_is_open(global::OFS::TRAJ)) {
            auto fname = utils::concat(save, "-", itraj, ".traj");
            utils::clearFile(fname);
            ofs_TRAJ.open(fname);
        }
    }
    return 0;
};

int Traj_Solver::init(const int& itraj) {
    init_ofs(itraj);

    if (itraj < 0) {  // pre-initialization
        pForceField->ForceField_init(nr0, np0, nm, N, itraj);
        for (int j = 0; j < N; ++j) nr[j] = 0.0f, np[j] = 0.0f;
    } else {  // initialization for each trajectory
        pForceField->ForceField_init(nr0, np0, nm, N, itraj);
        for (int j = 0; j < N; ++j) nr[j] = nr0[j], np[j] = np0[j];

        if (FLAGS_r != "") rst_read(itraj);
        LOG_IF(INFO, if_debug) << "initial traj done";
    }
    return 0;
};

int Traj_Solver::rst_read(const int& traj_in) {
    hload(pContext, "position", nr, N);
    hload(pContext, "momentum", np, N);
    return 0;
}

int Traj_Solver::rst_output(const int& traj_in) {
    hdump(pContext, "position", nr, N);
    hdump(pContext, "momentum", np, N);
    return 0;
}

int Traj_Solver::final(const int& itraj) {
    if (itraj < 0) {
        utils::closeOFS(ofs_SAMP);
    } else {
        utils::closeOFS(ofs_ENER);
        utils::closeOFS(ofs_TRAJ);
        rst_output(itraj);
    }
    return 0;
}

int Traj_Solver::check_break(int& succ) {
    if (succ != 0 || !ARRAY_ISFINITE(nr, N) || !ARRAY_ISFINITE(np, N)) succ = -1;
    return succ;
}

int Traj_Solver::traj_property(const num_real& dt) {
    Khere = 0.0f;
    for (int j = 0; j < N; ++j) Khere += 0.5f * np[j] * np[j] / nm[j];
    Vhere = vpes[0];
    Hhere = Khere + Vhere;
    Lhere = Khere - Vhere;
    Shere += Lhere * dt;
    return 0;
}

int Traj_Solver::traj(TCFnucl& tcfer, const int& N) { return traj_velocityverlet(tcfer, N); };


int Traj_Solver::update_p(const num_real& dt_in) {
    plFunction();
    for (int i = 0; i < N; ++i) np[i] -= nf[i] * dt_in;
    return 0;
}

int Traj_Solver::update_r(const num_real& dt_in) {
    plFunction();
    for (int i = 0; i < N; ++i) nr[i] += np[i] / nm[i] * dt_in;  //@slow to be opt (invM)
    return 0;
}

int Traj_Solver::update_thermo(const num_real& dt_in) {
    pThermo->evolve(nr, np, nm, dt_in, N);
    return 0;
}

int Traj_Solver::traj_velocityverlet(TCFnucl& tcfer, const int& N) {
    init(itraj);
    ff_calc1(level);
    LOG_IF(INFO, if_debug) << "cal force done";
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        LOG_IF(INFO, if_debug) << "before";
        if (succ == 0) succ = update_p(halfdt);
        LOG_IF(INFO, if_debug) << "p1";
        if (succ == 0) succ = update_r(halfdt);
        LOG_IF(INFO, if_debug) << "x1";
        if (succ == 0 && pThermo->dothermo(istep_dummy)) succ = update_thermo(dt);
        LOG_IF(INFO, if_debug) << "pthermo";
        if (succ == 0) update_r(halfdt);
        LOG_IF(INFO, if_debug) << "x2";
        if (succ == 0) succ = ff_calc1(level);
        LOG_IF(INFO, if_debug) << "f";
        // std::cout<< "3" << std::endl;
        if (succ == 0) update_p(halfdt);
        LOG_IF(INFO, if_debug) << "p2";
        if (succ == 0) traj_property(dt);

        if ((istep_dummy + 1) % sstep == 0) {
            if (check_break(succ) != 0) break;  // check only when do sampling
            // std::cout<< "4" << std::endl;
            // sampling
            isamp = (istep_dummy + 1) / sstep;
            if (succ == 0) succ = sampler(isamp, tcfer);
            // std::cout<< "5" << std::endl;
        }
    }
    final(itraj);
    return 0;
}

int Traj_Solver::sampler(const int& isamp, TCFnucl& tcfer) {
    estimator(isamp, tcfer);
    ispec = pForceField->ForceField_spec(nr, np, nm, N);
    // tcfer.Count(isamp, ispec);
    return 0;
}

int Traj_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
    num_real esti_V, esti_Kprim, esti_Kvir, esti_Eprim, esti_Evir, esti_Cprim, esti_Cvir;
    esti_V     = Vhere;
    esti_Kprim = 0.5f * N / pThermo->beta;
    esti_Kvir  = esti_Kprim;
    esti_Eprim = esti_V + esti_Kprim;
    esti_Evir  = esti_V + esti_Kvir;

    num_real sampunit = dt * sstep * iou.time;
    ofs_ENER << FMT(8) << itraj << FMT(8) << isamp * sampunit  //
             << FMT(8) << esti_V                               //
             << FMT(8) << esti_Kprim                           //
             << FMT(8) << esti_Kvir                            //
             << FMT(8) << esti_Eprim                           //
             << FMT(8) << esti_Evir                            //
             << std::endl;

    // ofs.open(std::to_string(itraj) + "trj.dat", std::ios_base::app);
    ofs_TRAJ << FMT(8) << isamp * sampunit;  //
    for (int i = 0; i < N; ++i) ofs_TRAJ << FMT(8) << nr[i];
    for (int i = 0; i < N; ++i) ofs_TRAJ << FMT(8) << np[i];
    for (int i = 0; i < N; ++i) ofs_TRAJ << FMT(8) << nf[i];
    ofs_TRAJ << std::endl;
    return 0;
}

int Traj_Solver::run_impl() {
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0

    TCFnucl coll;
    traj(coll, N);

    final(-1);
    return 0;
}

int Traj_Solver::run_parallel() {
    int nsave = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    if (nsave > FLAGS_nsave_mpi) nsave = FLAGS_nsave_mpi;
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0

    TCFnucl coll;

    for (int isave = 0; isave < nsave; ++isave) {
        int eachstart = (isave * ntraj) / nsave, eachend = ((isave + 1) * ntraj) / nsave, istart, iend;
        mpi_range(eachstart, eachend, mpi_nprocs, mpi_rank, istart, iend);
        CHECK_EQ(ntraj % (nsave * mpi_nprocs), 0);
        LOG(INFO) << "During [" << eachstart << ", " << eachend << "), "
                  << "mpi-" << mpi_rank << " cycle in [" << istart << ", " << iend << ")";

        MPI_Barrier(MPI_COMM_WORLD);
        for (int icycle = istart; icycle < iend; ++icycle) {
            itraj = icycle;
            // coll.Clear();
            try {
                traj(coll, N);
            } catch (std::runtime_error& e) {  // if some error, output currect results
                LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
            }
            // collsum.Amount(coll);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        // collmpi.MPIAmount(collsum);
        // if (mpi_rank == 0) collmpi.report(save + "-cache" + std::to_string(isave), sampunit);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // collmpi.MPIAmount(collsum);
    // if (mpi_rank == 0) collmpi.report(save, sampunit);
    // collsum.report(save + "-mpi" + std::to_string(mpi_rank), sampunit);

    final(-1);
    return 0;
}

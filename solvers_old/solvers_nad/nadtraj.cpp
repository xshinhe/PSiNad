#include "nadtraj.h"

#include <fstream>

#include "../solvers_md/traj.h"
#include "../utils/definitions.h"
#include "../utils/nad_utils.h"
#include "nadtcfer.h"

using namespace ARRAY_EG;

NadTraj_Solver::NadTraj_Solver(const Param& iparm, Model* pM) : Traj_Solver(iparm, pM) {
    // cast parameters
    pForceField = dynamic_cast<Nad_ForceField*>(pM);

    std::string mem_flag = Param_GetT(std::string, parm, "mem_flag", "#def");
    mem_type             = nad_mem::_dict.at(mem_flag);

    std::string eom_flag = Param_GetT(std::string, parm, "eom_flag", "#eac");
    eom_type             = elec_eom::_dict.at(eom_flag);

    std::string rep_flag = Param_GetT(std::string, parm, "rep_flag", "#dia");
    rep_type             = representation::_dict.at(rep_flag);

    std::string ini_flag = Param_GetT(std::string, parm, "ini_flag", "#occ");
    ini_type             = elec_init::_dict.at(ini_flag);

    std::string tcf_flag = Param_GetT(std::string, parm, "tcf_flag", "#rho");
    tcf_type             = nad_tcf::_dict.at(tcf_flag);

    Param_GetV(dyn_type, iparm, 0);  // -1, mdsampling; 0, real dynamics

    // get dimensions
    F    = pForceField->get_F();
    FF   = F * F;
    NF   = N * F;
    NFF  = N * F * F;
    NNF  = N * N * F;
    NNFF = N * N * F * F;
    FFFF = F * F * F * F;

    try {  // allocate memory
        // varibles for evolution of trajectory
        ///< time 0 & time t saves
        ALLOCATE_PTR_TO_VECTOR(eac0, F);
        ALLOCATE_PTR_TO_VECTOR(rho0, FF);
        ALLOCATE_PTR_TO_VECTOR(T0, FF);
        ALLOCATE_PTR_TO_VECTOR(rhot, FF);

        ///< dynamical variables for nuclear DOFs
        ALLOCATE_PTR_TO_VECTOR(fmean, N);
        ALLOCATE_PTR_TO_VECTOR(fcorr, N);
        ALLOCATE_PTR_TO_VECTOR(V, FF);
        ALLOCATE_PTR_TO_VECTOR(dV, NFF);
        ALLOCATE_PTR_TO_VECTOR(ddV, NNFF);

        // T: representation transform to adiabatic
        ALLOCATE_PTR_TO_VECTOR(T, FF);
        ALLOCATE_PTR_TO_VECTOR(E, F);
        ALLOCATE_PTR_TO_VECTOR(dE, NFF);
        ALLOCATE_PTR_TO_VECTOR(nacv, NFF);
        ALLOCATE_PTR_TO_VECTOR(ddE, NNFF);

        // S: representation transform to eigen-adiabatic
        ALLOCATE_PTR_TO_VECTOR(S, FF);
        ALLOCATE_PTR_TO_VECTOR(L, F);
        ALLOCATE_PTR_TO_VECTOR(dL, NFF);
        ALLOCATE_PTR_TO_VECTOR(ddL, NNFF);

        ///< dynamical variables for electronic DOFs
        ALLOCATE_PTR_TO_VECTOR(mvc, 2 * F);
        ALLOCATE_PTR_TO_VECTOR(eac, F);
        ALLOCATE_PTR_TO_VECTOR(eacf, F);
        ALLOCATE_PTR_TO_VECTOR(eacb, F);
        ALLOCATE_PTR_TO_VECTOR(rho, FF);
        ALLOCATE_PTR_TO_VECTOR(gmat, FF);
        ALLOCATE_PTR_TO_VECTOR(U, FF);
        ALLOCATE_PTR_TO_VECTOR(H, FF);
        ALLOCATE_PTR_TO_VECTOR(dH, NFF);
        ALLOCATE_PTR_TO_VECTOR(ddH, NNFF);

        // allocate temporary work array
        int max_work_size = std::max({2 * NFF, 4 * NF, 4 * FF, 2 * FF * FF});
        ALLOCATE_PTR_TO_VECTOR(workr, max_work_size);
        ALLOCATE_PTR_TO_VECTOR(workc, max_work_size);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }
}
NadTraj_Solver::~NadTraj_Solver(){};

int NadTraj_Solver::ff_calc1(const int& level,
                             const bool& refered) {  // forcefield calculation at a fixed nr, np
    plFunction();

    int succ = 0;
    switch (pForceField->type) {
        case ForceFieldModel:  // representation::diabatic & representation::adiabatic
            if (succ == 0) succ = pForceField->ForceField_npes(vpes, grad, hess, nr, np, level, N, itraj, istep);
            if (succ == 0) succ = pForceField->ForceField_epes(V, dV, ddV, nr, level, N, F, itraj, istep);
            break;
        case ForceFieldOnTheFly:
            if (rep_type != representation::onthefly) Param_Reset(rep_type, representation::onthefly);
            if (succ == 0) succ = pForceField->ForceField_epes(E, dE, ddE, nr, level, N, F, itraj, istep);
            break;
        default:
            LOG(FATAL) << "unknown ForceField type";
    }
    solve_transform(H, dH, ddH, S, L, dL, ddL, T, E, dE, ddE, V, dV, ddV, nr, np, nm, rep_type, level, N, F, workr,
                    workc, refered);
    return succ;
}

int NadTraj_Solver::ff_calc2() {
    plFunction();

    num_real* ForceMat;
    switch (rep_type) {
        case representation::diabatic:
            ForceMat = dV;  // using diabatic force matrix
            break;
        case representation::adiabatic:
        case representation::onthefly:
            ForceMat = dE;  // using adiabatic force matrix
            break;
        default:
            LOG(FATAL);
    }

    if (dyn_type == -1) {  // md sampling in ground state
        for (int j = 0; j < N; ++j) nf[j] = dE[j * FF];
        return 0;
    }

    switch (eom_type) {
        case elec_eom::eac:
        case elec_eom::eaccv:
            rho_eac(rho, eac, F);  // rho as temporary array
            for (int i = 0; i < FF; ++i) rho[i] -= gmat[i];
            break;
        case elec_eom::eacfb:
            rho_eac(rho, eacf, F);   // rho as temporary array
            rho_eac(gmat, eacb, F);  // gmat as temporary array
            for (int i = 0; i < FF; ++i) rho[i] = num_complex(0.5f) * (rho[i] + gmat[i]);
            break;
        case elec_eom::rho:  // gmat has been contained in rho
            break;
        default:
            LOG(FATAL);
    }
    {
        plScope("opt");
        // if you can do optimization to reduce the force, please override for the Nad_ForceField class
        pForceField->reduce_force(fmean, rho, ForceMat, N, F);
    }

    for (int i = 0; i < N; ++i) nf[i] = grad[i] + fmean[i];  // total force
    return 0;
}

int NadTraj_Solver::evolve_elec(num_complex* Uevolve) {
    plFunction();
    switch (eom_type) {
        case elec_eom::eac:
            update_eac(eac, Uevolve, N, F, workr, workc);
            break;
        case elec_eom::eaccv:
            update_eac(eac, Uevolve, N, F, workr, workc);
            update_rho(gmat, Uevolve, N, F, workr, workc);
            break;
        case elec_eom::eacfb:
            update_eac(eacf, Uevolve, N, F, workr, workc);
            update_eac(eacb, Uevolve, N, F, workr, workc);
            break;
        case elec_eom::rho:
            update_rho(rho, Uevolve, N, F, workr, workc);  // @NOTE: rho contains gmat!
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int NadTraj_Solver::init_occ2eac(const int& itraj) {
    for (int i = 0; i < F; ++i) mvc[i] = (i == occ0) ? 1 : 0;
    samp_mvc_focus(mvc, F);
    eac_mvc(eac0, mvc, F);
    return 0;
}

int NadTraj_Solver::init_ofs(const int& itraj) {
    if (itraj < 0) {
        if (global::ofs_is_open(global::OFS::SAMP)) {
            auto fname = utils::concat(save, "-", mpi_rank, ".samp");
            utils::clearFile(fname);
            ofs_SAMP.open(fname);
        }
        if (global::ofs_is_open(global::OFS::ESAMP)) {
            auto fname = utils::concat(save, "-", mpi_rank, ".esamp");
            utils::clearFile(fname);
            ofs_ESAMP.open(fname);
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
        if (global::ofs_is_open(global::OFS::ETRAJ)) {
            auto fname = utils::concat(save, "-", itraj, ".etraj");
            utils::clearFile(fname);
            ofs_ETRAJ.open(fname);
        }
    }
    return 0;
};

int NadTraj_Solver::init(const int& itraj) {
    init_ofs(itraj);

    plFunction();
    if (itraj < 0) {  // total inital
        // get information of: nm, occ0 and check
        pForceField->ForceField_init(nr0, np0, nm, rho0, eac0, occ0, N, F, itraj);
        for (int j = 0; j < N; ++j) nr[j] = 0.0f, np[j] = 0.0f;

        CHECK_GE(occ0, 0);
        CHECK_LT(occ0, F);

        // initialize T0
        switch (ini_type) {
            case elec_init::occ:
            case elec_init::eac:
            case elec_init::rho:
                ARRAY_EYE(T0, F);
                break;
            case elec_init::eig:
            case elec_init::d2a:
                ff_calc1(level);  // calculate forceField and solve eigen problem
                for (int i = 0; i < FF; ++i) T0[i] = T[i];
                break;
            case elec_init::rot: {
                CHECK_EQ(F, 2);
                num_real rotangle = Param_GetV(rotangle, parm, 0.0f) * phys::math::twopi;
                T0[0] = std::cos(rotangle), T0[1] = -std::sin(rotangle);
                T0[2] = std::sin(rotangle), T0[3] = std::cos(rotangle);
                break;
            }
            case elec_init::def: {
                try {
                    std::ifstream ifs("ini0");
                    num_real tmp;
                    for (int i = 0; i < FF; ++i)
                        if (ifs >> tmp) T0[i] = tmp;
                    ifs.close();
                } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }
                break;
            }
            default:
                LOG(FATAL);
        }
    } else {
        pForceField->ForceField_init(nr0, np0, nm, rho0, eac0, occ0, N, F, itraj);
        for (int j = 0; j < N; ++j) nr[j] = nr0[j], np[j] = np0[j];

        switch (ini_type) {
            case elec_init::occ:
            case elec_init::rot:
            case elec_init::def:
            case elec_init::eig:
            case elec_init::d2a:
                init_occ2eac(itraj);
                if (ini_type == elec_init::d2a) {  // update T0
                    ff_calc1(level);
                    for (int i = 0; i < FF; ++i) T0[i] = T[i];
                }
                ARRAY_MATMUL_TRANS1(eac, T0, eac0, F, F, 1);  // trans from dia T0-rep

                if (eom_type == elec_eom::eaccv) {  // @TODO: BUGS?
                    for (int i = 0, idx = 0; i < F; ++i)
                        for (int j = 0; j < F; ++j, ++idx)
                            gmat[idx] = (i == j) ? (NORM_OF(eac[i]) - ((occ0 == i) ? 1 : 0)) : phys::math::iz;
                    ARRAY_MATMUL_TRANS1(workc, T0, gmat, F, F, F);  // @fixed bug (Gadia = T^ Gdia T)
                    ARRAY_MATMUL(gmat, workc, T0, F, F, F);
                    // otherwise in later, gmat = gamma0*I, needn't to be transformed
                }
                break;
            case elec_init::eac: {
                try {
                    num_real tmp;
                    std::ifstream ifs("ini0");
                    for (int i = 0; i < F; ++i)
                        if (ifs >> tmp) eac[i] = tmp;
                    for (int i = 0; i < F; ++i)
                        if (ifs >> tmp) eac[i] += phys::math::im * tmp;
                    ifs.close();
                } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }
                break;
            }
            case elec_init::rho: {
                try {
                    num_real tmp;
                    std::ifstream ifs("ini0");
                    for (int i = 0; i < FF; ++i)
                        if (ifs >> tmp) rho[i] = tmp;
                    for (int i = 0; i < FF; ++i)
                        if (ifs >> tmp) rho[i] += phys::math::im * tmp;
                    ifs.close();
                    CHECK_EQ(eom_type, elec_eom::rho);
                } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }
                break;
            }
            default:
                LOG(FATAL);
        }
    }
    return 0;
};

int NadTraj_Solver::final(const int& itraj) {
    if (itraj < 0) {
        utils::closeOFS(ofs_SAMP);
        utils::closeOFS(ofs_ESAMP);
    } else {
        utils::closeOFS(ofs_ENER);
        utils::closeOFS(ofs_TRAJ);
        utils::closeOFS(ofs_ETRAJ);
    }
    return 0;
}

int NadTraj_Solver::rst_read(const int& traj_in) {
    hload(pContext, "position", nr, N);
    hload(pContext, "momentum", np, N);
    hload(pContext, "electron", eac, F);
    hload(pContext, "density", rho, FF);
    return 0;
}

int NadTraj_Solver::rst_output(const int& traj_in) {
    hdump(pContext, "position", nr, N);
    hdump(pContext, "momentum", np, N);
    hdump(pContext, "electron", eac, F);
    hdump(pContext, "density", rho, FF);
    return 0;
}

int NadTraj_Solver::check_break(int& succ) {
    plFunction();
    if (succ != 0 || !ARRAY_ISFINITE(nr, N) || !ARRAY_ISFINITE(np, N) || !ARRAY_ISFINITE(eac, F)) succ = -1;
    return succ;
}

int NadTraj_Solver::traj_property(const num_real& dt) {
    plFunction();
    Khere = 0.0f;
    for (int j = 0; j < N; ++j) Khere += 0.5f * np[j] * np[j] / nm[j];
    Vhere = vpes[0];
    Vhere += REAL_OF(ARRAY_TRACE2(rho, V, F, F));
    Hhere = Khere + Vhere;
    Lhere = Khere - Vhere;
    Shere += Lhere * dt;
    return 0;
}

int NadTraj_Solver::traj(NAD_TCFer& tcfer, const int& N, const int& F) {
    plFunction();
    return traj_velocityverlet(tcfer, N, F);
}

int NadTraj_Solver::traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) {
    plFunction();
    init(itraj);
    istep    = 0, ff_calc1(level), ff_calc2();
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        if (succ == 0) succ = update_p(halfdt);
        if (succ == 0) succ = update_r(halfdt);
        if (succ == 0 && dyn_type == -1 && pThermo->dothermo(istep_dummy)) succ = update_thermo(dt);
        if (succ == 0) succ = update_r(halfdt);
        if (succ == 0) succ = ff_calc1(level, true);
        if (succ == 0) succ = solve_Ut(U, S, L, T, E, dt, rep_type, N, F, workr, workc);
        if (succ == 0) succ = evolve_elec(U);
        if (succ == 0) succ = ff_calc2();
        if (succ == 0) succ = update_p(halfdt);
        if (succ == 0) traj_property(dt);

        if ((istep_dummy + 1) % sstep == 0) {
            plScope("sample");
            isamp = (istep_dummy + 1) / sstep;
            if (check_break(succ) != 0) break;
            if (succ == 0) succ = sampler(isamp, tcfer);
            if (global::ofs_is_open(global::OFS::TRAJ)) {
                if (global::ofs_is_open(global::OFS::ETRAJ)) {
                    pForceField->ForceField_write(ofs_TRAJ, ofs_ETRAJ,  //
                                                  nr, np, nm, rhot, eac, occt, N, F, itraj, istep);
                } else {
                    pForceField->ForceField_write(ofs_TRAJ, ofs_TRAJ,  //
                                                  nr, np, nm, rhot, eac, occt, N, F, itraj, istep);
                }
            }
        }
    }

    final(itraj);
    return 0;
}

int NadTraj_Solver::sampler(const int& isamp, NAD_TCFer& tcfer) {
    plScope("sampler");
    correlation(isamp, tcfer);
    ispec = pForceField->ForceField_spec(nr, np, nm, N, F);
    tcfer.Count(isamp, ispec);
    return 0;
}

int NadTraj_Solver::kernel0(num_complex* rhox, const int& F) {
    LOG(FATAL);
    return 0;
}
int NadTraj_Solver::kernelt(num_complex* rhox, const int& F) {
    LOG(FATAL);
    return 0;
}
int NadTraj_Solver::correlation(const int& isamp, NAD_TCFer& tcfer) {
    if (isamp == 0) kernel0(rho0, F);
    kernelt(rhot, F);
    if (tcfer.tcf_reduced) {
        tcfer.val[0] = phys::math::iz;
        for (int rcdidx = 0, i0 = 0; i0 < FF; ++i0)
            for (int it = 0; it < FF; ++it, ++rcdidx)
                if (tcfer.tcf_0t_bool[rcdidx]) {
                    tcfer.val[0] += tcf_weight * tcfer.tcf_0t_val[rcdidx] * rho0[i0] * rhot[it];
                }
    } else {
        plScope("tcf");
        for (int idx = 0, rcdidx = 0, i0 = 0; i0 < FF; ++i0)
            for (int it = 0; it < FF; ++it, ++rcdidx)
                if (tcfer.tcf_0t_bool[rcdidx]) tcfer.val[idx++] = tcf_weight * rho0[i0] * rhot[it];
    }
    return 0;
}

int NadTraj_Solver::run_impl() {
    plFunction();

    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0
    NAD_TCFer coll = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);

    try {
        traj(coll, N, F);
    } catch (std::runtime_error& e) {  // if some error, output currect results
        LOG(WARNING) << "runtime_error cause breakdown for traj";
    }
    coll.report(save, sampunit);
    final(-1);
    return 0;
}


int NadTraj_Solver::run_parallel() {
    plFunction();

    int nsave = (ntraj > mpi_nprocs) ? ntraj / mpi_nprocs : 1;
    if (nsave > FLAGS_nsave_mpi) nsave = FLAGS_nsave_mpi;
    num_real sampunit = dt * sstep * iou.time;
    init(-1);  // get occ0
    NAD_TCFer coll    = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collsum = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);
    NAD_TCFer collmpi = NAD_TCFer(tcf_type, nspec, nsamp, N, F, occ0);

    for (int isave = 0; isave < nsave; ++isave) {
        int eachstart = (isave * ntraj) / nsave, eachend = ((isave + 1) * ntraj) / nsave, istart, iend;
        mpi_range(eachstart, eachend, mpi_nprocs, mpi_rank, istart, iend);
        CHECK_EQ(ntraj % (nsave * mpi_nprocs), 0);
        LOG(INFO) << "During [" << eachstart << ", " << eachend << "), "
                  << "mpi-" << mpi_rank << " cycle in [" << istart << ", " << iend << ")";

        MPI_Barrier(MPI_COMM_WORLD);
        for (int icycle = istart; icycle < iend; ++icycle) {
            itraj = icycle;
            coll.Clear();
            try {
                traj(coll, N, F);
            } catch (std::runtime_error& e) {  // if some error, output currect results
                LOG(WARNING) << "runtime_error cause breakdown for traj=" << itraj;
                // continue; // TODO ?? FIXBUG
            }
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
}

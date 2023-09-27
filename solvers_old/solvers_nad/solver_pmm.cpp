#include "solver_pmm.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

/**
 * @brief Perturbation Meyer-Miller Solver
 * @details the phase angle is integrated before
 * @param iparm [parameter structure]
 * @param pM [Model class]
 * @param M [description]
 * @param t [description]
 */
PMM_Solver::PMM_Solver(const Param& iparm, Model* pM)
    : Multi_NadTraj_Solver(iparm, pM, 2 * dynamic_cast<Nad_ForceField*>(pM)->get_F()) {
    /// pass it as M = 2*F

    restep = Param_GetT(int, parm, "restep", 20);

    Param_GetV(scale, parm, 1.0f);
    eom_type = elec_eom::rho;  // use rho for EOM

    try {
        ALLOCATE_PTR_TO_VECTOR(drhos, MFF);
        ALLOCATE_PTR_TO_VECTOR(deltarhos, MFF);
        ALLOCATE_PTR_TO_VECTOR(rho_corr, FF);
        ALLOCATE_PTR_TO_VECTOR(rhosumt, FF);
        ALLOCATE_PTR_TO_VECTOR(rhoavg, FF);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    save = name() + "_" + pForceField->tag;
}

PMM_Solver::~PMM_Solver(){};

int PMM_Solver::traj(NAD_TCFer& tcfer, const int& N, const int& F) {
    plFunction();
    init(itraj);
    ff_calc1(false, 1), ff_calc2();
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        if (succ == 0) succ = update_p(halfdt);
        if (succ == 0) succ = update_r(halfdt);
        if (succ == 0 && dyn_type == -1 && pThermo->dothermo(istep_dummy)) succ = update_thermo(dt);
        if (succ == 0) succ = update_r(halfdt);
        if (succ == 0) succ = ff_calc1(true, 1);
        // if (succ == 0) succ = solve_Ut(U, S, L, T, E, dt, rep_type, N, F, workr, workc);
        if (succ == 0) succ = evolve_elec(Us);
        if (succ == 0) succ = multi_ff_calc2();
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

int PMM_Solver::multi_evolve_elec(num_complex* Uevolve) {
    plFunction();

    num_complex* rho_0      = rhos;
    num_complex* rho_2      = rhos + (2 * F - 1) * FF;
    num_complex* drho_0     = drhos;
    num_complex* drho_2     = drhos + (2 * F - 1) * FF;
    num_complex* rhotmp     = deltarhos;
    num_complex* deltarho_2 = deltarhos + (2 * F - 1) * FF;
    num_complex* U_0        = Us;
    num_complex* U_2        = Us + (2 * F - 1) * FF;


    isumt += 1;
    for (int i = 0; i < FF; ++i) rhosumt[i] += rho_2[i] * dt, rhoavg[i] = rhosumt[i] / (isumt * dt);

    // drho0: Ehrenfest Path
    update_drho(drho_0, rho_0, U_0, N, F, workr, workc);

    // drho2: first step
    update_drho(drho_2, rho_0, U_2, N, F, workr, workc);
    // {
    //     std::ofstream ofs("rho2-a", std::ios_base::app);
    //     ofs << FMT(8) << istep * dt;
    //     for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(drho_2[i]) / dt;
    //     ofs << std::endl;
    //     ofs.close();
    // }

    update_drho(rhotmp, deltarho_2, U_0, N, F, workr, workc);
    // if (REAL_OF(rho_2[0]) < 1 && REAL_OF(rho_2[0]) > 0) {
    for (int i = 0; i < FF; ++i) drho_2[i] += rhotmp[i];  // - 10.0e0 * (rho_2[i] - rhoavg[i]) * dt;
    // }
    // {
    //     std::ofstream ofs("rho2-b", std::ios_base::app);
    //     ofs << FMT(8) << istep * dt;
    //     for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(drho_2[i]) / dt;
    //     ofs << std::endl;
    //     ofs.close();
    // }

    // drho1
    for (int i = 1; i < 2 * F - 1; ++i) {
        num_complex* rho_1      = rhos + i * FF;
        num_complex* drho_1     = drhos + i * FF;
        num_complex* deltarho_1 = deltarhos + i * FF;
        num_complex* U_1        = Us + i * FF;

        update_drho(drho_1, rho_0, U_1, N, F, workr, workc);
        update_drho(rhotmp, deltarho_1, U_0, N, F, workr, workc);
        for (int i = 0; i < FF; ++i) drho_1[i] += rhotmp[i];

        // drho2: second step
        update_drho(rhotmp, rho_1, U_1, N, F, workr, workc);
        for (int i = 0; i < FF; ++i) drho_2[i] += rhotmp[i] - drho_1[i];
        // {
        //     std::ofstream ofs(utils::concat("rho1-", i), std::ios_base::app);
        //     ofs << FMT(8) << istep * dt;
        //     for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(drho_1[i]) / dt;
        //     ofs << std::endl;
        //     ofs.close();
        // }
        // {
        //     std::ofstream ofs(utils::concat("rho2-c", i), std::ios_base::app);
        //     ofs << FMT(8) << istep * dt;
        //     for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(drho_2[i]) / dt;
        //     ofs << std::endl;
        //     ofs.close();
        // }
    }

    // collect & update
    for (int i = 0; i < M * FF; ++i) rhos[i] += drhos[i];
    for (int mu = 0, idx = 0; mu < M; ++mu) {
        for (int i = 0; i < FF; ++i, ++idx) deltarhos[idx] = rhos[idx] - rhos[i];
    }

    // {
    //     std::ofstream ofs("rhos", std::ios_base::app);
    //     ofs << FMT(4) << istep * dt;
    //     for (int i = 0; i < M * FF; ++i) ofs << FMT(8) << REAL_OF(rhos[i]);
    //     ofs << "\n";
    //     ofs.close();
    // }


    // collapse to ehrenfest again if necessary

    // if ((istep + 1) % restep == 0) multi_reinit(itraj);

    // std::ofstream ofs("test_r", std::ios_base::app);
    // for (int mu = 0, idx = 0; mu < M; ++mu) {
    //     for (int j = 0; j < N; ++j, ++idx) ofs << FMT(8) << nrs[idx];
    //     ofs << std::endl;
    // }
    // ofs.close();

    // ofs.open("test_p", std::ios_base::app);
    // for (int mu = 0, idx = 0; mu < M; ++mu) {
    //     for (int j = 0; j < N; ++j, ++idx) ofs << FMT(8) << nps[idx];
    //     ofs << std::endl;
    // }
    // ofs.close();

    // ofs.open("test_rho", std::ios_base::app);
    // for (int mu = 0, idx = 0; mu < M; ++mu) {
    //     for (int i = 0; i < FF; ++i, ++idx) ofs << FMT(8) << REAL_OF(rhos[idx]);
    //     ofs << std::endl;
    // }
    // ofs.close();

    return 0;
}

int PMM_Solver::multi_ff_calc2() {
    num_real* ForceMats  = (rep_type == representation::diabatic) ? dVs : dEs;
    num_real* ForceMat_0 = ForceMats;
    num_real* nf_0       = nfs;
    num_real* grad_0     = grads;
    num_complex* rho_0   = rhos;

    pForceField->reduce_force(fmean, rho_0, ForceMat_0, N, F);
    for (int i = 0; i < N; ++i) nf_0[i] = grad_0[i] + fmean[i];  // total force

    // calculate correllated force
    for (int mu = 1; mu < M; ++mu) {
        num_real* ForceMat_mu    = ForceMats + mu * NFF;
        num_real* nf_mu          = nfs + mu * N;
        num_real* grad_mu        = grads + mu * N;
        num_complex* deltarho_mu = deltarhos + mu * FF;
        for (int i = 0; i < N; ++i) nf_mu[i] = grad_mu[i];
        pForceField->reduce_force(fmean, rho_0, ForceMat_mu, N, F);
        for (int i = 0; i < N; ++i) nf_mu[i] += fmean[i];
        pForceField->reduce_force(fmean, deltarho_mu, ForceMat_0, N, F);
        for (int i = 0; i < N; ++i) nf_mu[i] += fmean[i];
    }
    return 0;
}

int PMM_Solver::multi_init(const int& itraj) {
    NadTraj_Solver::init(itraj);

    // for (int i = 0; i < N; ++i) nr[i] = 0.01e0, np[i] = -0.01e0;

    if (ini_type == elec_init::occ || ini_type == elec_init::d2a) {
        // copy nuclear variables
        for (int mu = 0, idx = 0; mu < M; ++mu) {
            for (int j = 0; j < N; ++j, ++idx) nrs[idx] = nr[j], nps[idx] = np[j], nms[idx] = nm[j];
        }

        rho_eac(rho, eac0, F);
        for (int mu = 0, idx = 0; mu < M; ++mu) {
            for (int i = 0; i < FF; ++i, ++idx) rhos[idx] = rho[i];
        }

        // add perturbations to rho(1)
        for (int k = 0, idx = FF; k < F; ++k) {
            if (k == occ0) continue;
            num_complex* rhox = rhos + idx;
            rhox[k * F + occ0] += 0.5e0 * scale;
            rhox[occ0 * F + k] += 0.5e0 * scale;
            idx += FF;
            num_complex* rhoy = rhos + idx;
            rhoy[k * F + occ0] += 0.5e0 * phys::math::im * scale;
            rhoy[occ0 * F + k] -= 0.5e0 * phys::math::im * scale;
            idx += FF;
        }

        if (ini_type == elec_init::d2a) {
            ff_calc1(false, 1);
            for (int mu = 0; mu < M; ++mu) {
                num_real* T_mu       = Ts + mu * FF;
                num_complex* rhos_mu = rhos + mu * FF;
                ARRAY_MATMUL_TRANS1(workc, T_mu, rhos_mu, F, F, F);
                ARRAY_MATMUL(rhos_mu, workc, T_mu, F, F, F);
            }
        }

        // init deltarhos
        for (int mu = 0, idx = 0; mu < M; ++mu) {
            for (int i = 0; i < FF; ++i, ++idx) deltarhos[idx] = rhos[idx] - rhos[i];
        }
        // clean rho_corr
        memset(rho_corr, 0, FF * sizeof(num_complex));
        memset(rhosumt, 0, FF * sizeof(num_complex));
        isumt = 0;
    } else {
        LOG(FATAL);
    }

    return 0;
};


int PMM_Solver::multi_reinit(const int& itraj) {
    // ARRAY_SHOW(rhos, M, FF);

    // backup previsous perturbations
    // num_complex* rho_2 = rhos + (2 * F - 1) * FF;
    // for (int i = 0; i < FF; ++i) rho_corr[i] += (rho_2[i] - rhos[i]) / (scale * scale);

    // copy Ehrenfest's nuclear/electronic variables
    for (int mu = 0, idx = 0; mu < M; ++mu) {
        for (int j = 0; j < N; ++j, ++idx) nrs[idx] = nrs[j], nps[idx] = nps[j], nms[idx] = nms[j];
    }

    return 0;

    for (int mu = 0, idx = 0; mu < M - 1; ++mu) {
        for (int i = 0; i < FF; ++i, ++idx) rhos[idx] = rhos[i];
    }

    num_complex* Teac = workc;
    EigenSolve(workr, Teac, rhos, F);  // solve eigen eac for ehrenfest's rho

    // ARRAY_SHOW(workr, 1, F);
    // ARRAY_SHOW(Teac, F, F);

    // add perturbations to rho(1)
    for (int k = 0, idx = FF, km = 0; k < F; ++k) {
        if (k == occ0) continue;
        num_complex* rhox = rhos + idx;

        for (int i1 = 0, i12 = 0; i1 < F; ++i1) {
            for (int i2 = 0; i2 < F; ++i2, ++i12) {
                rhox[i12] += 0.5e0 * scale *
                             (Teac[i1 * F + km] * Teac[i2 * F + F - 1] + Teac[i2 * F + km] * Teac[i1 * F + F - 1]);
            }
        }
        idx += FF;

        num_complex* rhoy = rhos + idx;
        for (int i1 = 0, i12 = 0; i1 < F; ++i1) {
            for (int i2 = 0; i2 < F; ++i2, ++i12) {
                rhoy[i12] += 0.5e0 * scale * phys::math::im *
                             (Teac[i1 * F + km] * Teac[i2 * F + F - 1] - Teac[i2 * F + km] * Teac[i1 * F + F - 1]);
            }
        }
        idx += FF;
        km++;
    }
    // init deltarhos
    for (int mu = 0, idx = 0; mu < M - 1; ++mu) {
        for (int i = 0; i < FF; ++i, ++idx) deltarhos[idx] = rhos[idx] - rhos[i];
    }

    // ARRAY_SHOW(rhos, M, FF);
    // LOG(FATAL);

    return 0;
};

int PMM_Solver::kernel0(num_complex* rhox, const int& F) {
    num_complex* rho_0 = rhos;
    num_complex* rho_2 = rhos + (2 * F - 1) * FF;

    for (int i = 0; i < FF; ++i) rhox[i] = rho_0[i] + 1.0e0 * (rho_2[i] - rho_0[i]) / (scale * scale);

    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        num_real* T_0 = Ts;
        num_real* T_2 = Ts + (2 * F - 1) * FF;
        ARRAY_MATMUL3_TRANS2(rhox, T_0, rho_0, T_0, F, F, F, F);
        ARRAY_MATMUL3_TRANS2(workc, T_2, rho_2, T_2, F, F, F, F);
        // for (int i = 0; i < FF; ++i) rhox[i] = (1 - 1.0e0 / (scale * scale)) * rhox[i] + workc[i] / (scale * scale);

        // ARRAY_MATMUL3_TRANS2(workc, T_0, rhox, T_0, F, F, F, F);
        // for (int i = 0; i < FF; ++i) rhox[i] = workc[i];
    }
    return 0;
}

int PMM_Solver::kernelt(num_complex* rhox, const int& F) { return kernel0(rhox, F); }

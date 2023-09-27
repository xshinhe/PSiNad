#include "solver_wmm.h"

#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

WMM_Solver::WMM_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    xi0 = 1.0f;
    Param_GetV(gamma0, iparm, 0.05);
    Param_GetV(adjustp, iparm, 0);
    gammat     = gamma0;
    totact     = (1 + F * gamma0) / xi0;
    xit        = xi0;
    tcf_weight = num_complex(F);

    std::string suffix;
    suffix += ((rep_type == representation::adiabatic) ? "adia" : "");
    suffix += ((eom_type == elec_eom::eaccv) ? "cv" : "");
    suffix += std::to_string(gamma0);
    save = name() + suffix + "_" + pForceField->tag;
};

WMM_Solver::~WMM_Solver(){};

int WMM_Solver::init_occ2eac(const int& itraj) {
    samp_mvc_sphere(mvc, 2.0f * totact, F);
    eac_mvc(eac0, mvc, F);
    return 0;
}

int WMM_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);
    if (eom_type != elec_eom::eaccv) {  // override gmat
        for (int i = 0, idx = 0; i < F; ++i)
            for (int j = 0; j < F; ++j, ++idx) gmat[idx] = (i == j) ? gamma0 : phys::math::iz;
    }
    // calc kernel
    if (adjustp == 1) {
        ff_calc1(level);
        kernel0(rho0, F);
        num_real Ekin = 0.0f;
        num_real Epot = REAL_OF(ARRAY_TRACE2(rho0, V, F, F));
        num_real Eidl = V[occ0 * (F + 1)];
        for (int i = 0; i < N; ++i) { Ekin += 0.5f * np[i] * np[i] / nm[i]; }
        num_real newEkin = Ekin + Eidl - Epot;
        for (int i = 0; i < N; ++i) { np[i] *= sqrt(newEkin / Ekin); }
    }

    return 0;
}
int WMM_Solver::kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F) {
    rho_eac(rhox, eac, F);
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    for (int idx = 0, i = 0; i < F; ++i)
        for (int j = 0; j < F; ++j, ++idx)
            rhox[idx] = xic * rhox[idx] - ((i == j) ? gammac * phys::math::iu : phys::math::iz);
    return 0;
}

int WMM_Solver::kernel0(num_complex* rhox, const int& F) { return kernel_cmm(rhox, xi0, gamma0, F); }
int WMM_Solver::kernelt(num_complex* rhox, const int& F) { return kernel_cmm(rhox, xit, gammat, F); }

int WMM_Solver::traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) {
    init(itraj);
    ff_calc1(level), ff_calc2();
    int succ = sampler(0, tcfer);
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep
        if (succ == 0)
            for (int i = 0; i < N; ++i) np[i] -= nf[i] * halfdt;
        if (succ == 0)
            for (int i = 0; i < N; ++i) nr[i] += np[i] / nm[i] * halfdt;
        if (succ == 0 && dyn_type == -1 && pThermo->dothermo(istep_dummy)) succ = pThermo->evolve(nr, np, nm, dt, N);
        if (succ == 0)
            for (int i = 0; i < N; ++i) nr[i] += np[i] / nm[i] * halfdt;
        if (succ == 0) succ = ff_calc1(level);
        if (succ == 0) succ = solve_Ut(U, S, L, T, E, dt, rep_type, N, F, workr, workc);
        if (succ == 0) succ = evolve_elec(U);
        if (succ == 0) succ = ff_calc2();
        if (succ == 0)
            for (int i = 0; i < N; ++i) np[i] -= nf[i] * halfdt;
        if (succ == 0) traj_property(dt);

        if ((istep_dummy + 1) % sstep == 0) {
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

    // LOG(WARNING);

    if (ofs_SAMP.is_open()) {
        ofs_SAMP << FMT(0) << 0 << FMT(0) << itraj;
        ofs_SAMP << FMT(8) << (istep + 1) * dt  //
                 << FMT(8) << nr[0]             //
                 << FMT(8) << np[0]             //
                 << FMT(8) << nr0[0]            //
                 << FMT(8) << np0[0];
        for (int i = 0; i < tcfer.lentcf; ++i)
            ofs_SAMP << FMT(8) << REAL_OF(tcfer.val[i]) << FMT(8) << IMAG_OF(tcfer.val[i]);
        for (int i = 0; i < FF; ++i) ofs_SAMP << FMT(8) << F * REAL_OF(rho0[i]) << FMT(8) << F * IMAG_OF(rho0[i]);
        ofs_SAMP << std::endl;
    } else if (itraj == 0) {
        LOG(WARNING) << "NOT OPEN";
    }

    final(itraj);
    return 0;
}

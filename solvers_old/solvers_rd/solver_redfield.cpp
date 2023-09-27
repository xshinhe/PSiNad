#include "solver_redfield.h"

#include "../utils/definitions.h"

using namespace ARRAY_EG;

RedField_Solver::RedField_Solver(Param iparm, Model* pM) : Solver(iparm, pM) {
    pForceField = dynamic_cast<SystemBath_ForceField*>(pM);

    F     = pForceField->get_F();
    N     = pForceField->get_N();
    nbath = pForceField->get_nbath();
    Nb    = pForceField->get_Nb();

    FF = F * F, NN = N * N, NF = N * F;
    NFF = N * F * F, NNF = N * N * F;
    NNFF = N * N * F * F, FFFF = F * F * F * F;
    NbFF = Nb * FF;

    Param_GetV(tend, iparm, 0.0f);
    Param_GetV(dt, iparm, -1.0f);
    Param_GetV(sstep, iparm, 1);
    if (pForceField->Suggest_tend() > 0.0f) tend = pForceField->Suggest_tend();
    if (pForceField->Suggest_dt() > 0.0f) dt = pForceField->Suggest_dt();

    vscale = 1;
    if (dt > 0.1f) {
        LOG(WARNING) << "Large dt will cause unconvergece for RK-4 expansion, so do scaling on dt=0.1";
        vscale = 0.01f / dt;   // remembering [M] ~ [L-1] ~ [T-1]
        dt     = vscale * dt;  // = 0.01
        tend   = vscale * tend;
    }
    nstep = sstep * (int(tend / (sstep * dt)) + 1);  // make to over off
    nsamp = nstep / sstep + 1;
    CHECK_GT(tend, 0);
    CHECK_GT(dt, 0);
    CHECK_GT(sstep, 0);
    CHECK_GT(nstep, 0);
    CHECK_GT(nsamp, 0);

    try {
        ALLOCATE_PTR_TO_VECTOR(Eele, F);
        ALLOCATE_PTR_TO_VECTOR(Tele, FF);
        ALLOCATE_PTR_TO_VECTOR(Hele, FF);
        ALLOCATE_PTR_TO_VECTOR(mDE, FF);
        ALLOCATE_PTR_TO_VECTOR(C_mDE, FF);

        ALLOCATE_PTR_TO_VECTOR(Qtran, nbath * FF);
        ALLOCATE_PTR_TO_VECTOR(GM_tensor, FF * FF);
        ALLOCATE_PTR_TO_VECTOR(R_tensor, FF * FF);

        ALLOCATE_PTR_TO_VECTOR(eac0, F);
        ALLOCATE_PTR_TO_VECTOR(rho0, FF);
        ALLOCATE_PTR_TO_VECTOR(rho1, FF);
        ALLOCATE_PTR_TO_VECTOR(rho2, FF);
        ALLOCATE_PTR_TO_VECTOR(rho3, FF);
        ALLOCATE_PTR_TO_VECTOR(rho4, FF);
        ALLOCATE_PTR_TO_VECTOR(rhodia, FF);
        ALLOCATE_PTR_TO_VECTOR(rhoadia, FF);
        ALLOCATE_PTR_TO_VECTOR(rhotmp, FF);

        ALLOCATE_PTR_TO_VECTOR(workr, FF);
        ALLOCATE_PTR_TO_VECTOR(workc, FF);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }

    // calc adiabatic basis
    for (int i = 0; i < FF; ++i) Hele[i] = pForceField->Hsys[i];
    EigenSolve(Eele, Tele, Hele, F);  // solve and save pure electronic probelm
    for (int i = 0, ik = 0; i < F; ++i) {
        for (int k = 0; k < F; ++k, ++ik) mDE[ik] = -(Eele[i] - Eele[k]);
    }

    double* omegas = pForceField->omegas;  // temporary
    double* coeffs = pForceField->coeffs;  // temporary

    pForceField->mybath->fun_Cw(C_mDE, mDE, FF, omegas, coeffs, Nb);

    // calc Qtran in adiabatic basis
    for (int ib = 0; ib < nbath; ++ib) {
        double* Qtrani = Qtran + ib * FF;
        double* Qi     = pForceField->Q + ib * FF;
        ARRAY_MATMUL3_TRANS1(Qtrani, Tele, Qi, Tele, F, F, F, F);
    }
    calc_R_tensor();

    // so scaling for energy
    for (int i = 0; i < FF; ++i) mDE[i] /= vscale, Hele[i] /= vscale;
    for (int i = 0; i < FFFF; ++i) GM_tensor[i] /= vscale, R_tensor[i] /= vscale;

    save = name() + "_" + pForceField->tag;
}

RedField_Solver::~RedField_Solver(){};

int RedField_Solver::calc_R_tensor() {
    // first calc GM_tensor
    for (int n1 = 0, n1m1 = 0, n1m1n2m2 = 0; n1 < F; ++n1) {
        for (int m1 = 0; m1 < F; ++m1, ++n1m1) {
            for (int n2 = 0, n2m2 = 0; n2 < F; ++n2) {
                for (int m2 = 0; m2 < F; ++m2, ++n2m2, ++n1m1n2m2) {
                    GM_tensor[n1m1n2m2] = 0.0e0;
                    for (int ib = 0, idxQi = 0; ib < nbath; ++ib, idxQi += FF) {
                        double* Qi = Qtran + idxQi;
                        for (int jb = 0, idxQj = 0; jb < nbath; ++jb, idxQj += FF) {
                            double* Qj = Qtran + idxQj;
                            if (ib == jb) {  // only when ib = jb, C_mDE != 0
                                GM_tensor[n1m1n2m2] += Qi[n1m1] * Qj[n2m2] * C_mDE[n2m2];
                            }
                        }
                    }
                }
            }
        }
    }
    // second calc Rtensor
    int FFF = FF * F;
    for (int n1 = 0, n1m1 = 0, n1m1n2m2 = 0; n1 < F; ++n1) {
        for (int m1 = 0; m1 < F; ++m1, ++n1m1) {
            for (int n2 = 0, n2m2 = 0; n2 < F; ++n2) {
                for (int m2 = 0; m2 < F; ++m2, ++n2m2, ++n1m1n2m2) {
                    R_tensor[n1m1n2m2] = GM_tensor[m2 * FFF + m1 * FF + n1 * F + n2] +
                                         CONJ_OF(GM_tensor[n2 * FFF + n1 * FF + m1 * F + m2]);
                    if (m1 == m2) {
                        for (int k = 0; k < F; ++k) { R_tensor[n1m1n2m2] -= GM_tensor[n1 * FFF + k * FF + k * F + n2]; }
                    }
                    if (n1 == n2) {
                        for (int k = 0; k < F; ++k) {
                            R_tensor[n1m1n2m2] -= CONJ_OF(GM_tensor[m1 * FFF + k * FF + k * F + m2]);
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int RedField_Solver::action_on_wavafunction(num_complex* rhonew, num_complex* rhoold, const num_real& t) {
    plFunction();
    ARRAY_MATMUL(rhonew, R_tensor, rhoold, FF, FF, 1);
    for (int i = 0; i < FF; ++i) rhonew[i] += phys::math::im * mDE[i] * rhoold[i];
    return 0;
}

int RedField_Solver::run_impl() {
    plFunction();

    // initialization
    pForceField->ForceField_init_elec(rho0, eac0, occ0, F, itraj);
    if (occ0 >= 0 && occ0 < F) {
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) rho0[ik] = (i == k) && (i == occ0) ? 1.0e0 : 0.0e0;
        }
    }
    ARRAY_MATMUL3_TRANS1(rhoadia, Tele, rho0, Tele, F, F, F, F);

    std::ofstream ofs(save);
    num_real timeunit = iou.time / vscale;
    num_real h1 = dt / 6.0f, h2 = dt / 3.0f, h3 = dt / 3.0f, h4 = dt / 6.0f;

    // time = 0
    ARRAY_MATMUL3_TRANS2(rhodia, Tele, rhoadia, Tele, F, F, F, F);

    ofs << FMT(8) << "stat" << FMT(8) << "time";
    for (int i1 = 0; i1 < F; ++i1)
        for (int i2 = 0; i2 < F; ++i2) {
            ofs << FMT(8) << utils::concat("Re(", i1, ",", i2, ")");
            ofs << FMT(8) << utils::concat("Im(", i1, ",", i2, ")");
        }
    ofs << std::endl;

    ofs << FMT(8) << 0 << FMT(8) << 0.0e0;
    for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(rhodia[i]) << FMT(8) << IMAG_OF(rhodia[i]);
    ofs << std::endl;

    int succ = 0;
    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;
        {
            if (succ == 0) succ = action_on_wavafunction(rho1, rhoadia, istep * dt);
            for (int i = 0; i < FF; ++i) rhotmp[i] = rhoadia[i] + 0.5f * dt * rho1[i];
            if (succ == 0) succ = action_on_wavafunction(rho2, rhotmp, (istep + 0.5f) * dt);
            for (int i = 0; i < FF; ++i) rhotmp[i] = rhoadia[i] + 0.5f * dt * rho2[i];
            if (succ == 0) succ = action_on_wavafunction(rho3, rhotmp, (istep + 0.5f) * dt);
            for (int i = 0; i < FF; ++i) rhotmp[i] = rhoadia[i] + 1.0f * dt * rho3[i];
            if (succ == 0) succ = action_on_wavafunction(rho4, rhotmp, (istep + 1.0f) * dt);
            for (int i = 0; i < FF; ++i) rhoadia[i] += (rho1[i] * h1 + rho2[i] * h2 + rho3[i] * h3 + rho4[i] * h4);
        }
        if ((istep_dummy + 1) % sstep == 0) {
            ARRAY_MATMUL3_TRANS2(rhodia, Tele, rhoadia, Tele, F, F, F, F);
            ofs << FMT(8) << 0 << FMT(8) << istep * dt * timeunit;
            for (int i = 0; i < FF; ++i) ofs << FMT(8) << REAL_OF(rhodia[i]) << FMT(8) << IMAG_OF(rhodia[i]);
            ofs << std::endl;
        }
    }
    return 0;
}

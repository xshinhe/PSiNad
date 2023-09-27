#include "nad_utils.h"

#include "definitions.h"

using namespace ARRAY_EG;

int eac_mvc(num_complex* eac, num_real* mvc, int fdim) {
    plFunction();
    num_real *mapvarx = mvc, *mapvarp = mvc + fdim;
    for (int i = 0; i < fdim; ++i) { eac[i] = phys::math::sqrthalf * (mapvarx[i] + phys::math::im * mapvarp[i]); }
    return 0;
}

int mvc_eac(num_real* mvc, num_complex* eac, int fdim) {
    plFunction();
    num_real *mapvarx = mvc, *mapvarp = mvc + fdim;
    for (int i = 0; i < fdim; ++i) {
        mapvarx[i] = phys::math::sqrttwo * REAL_OF(eac[i]);
        mapvarp[i] = phys::math::sqrttwo * IMAG_OF(eac[i]);
    }
    return 0;
}

int rho_eac(num_complex* rho, num_complex* eac, int fdim) {
    plFunction();
    ARRAY_OUTER_TRANS2(rho, eac, eac, fdim, fdim);
    return 0;
}

int samp_mvc_focus(num_real* mvc, int fdim) {  // assume that left half mvc saved actions
    plFunction();
    num_real* mvc_mid = mvc + fdim;
    rand_uniform(mvc_mid, fdim, phys::math::twopi);  // use right half mvc as random angles
    for (int i = 0; i < fdim; ++i) {
        mvc[i]     = phys::math::sqrttwo * sqrt(mvc[i]) * cos(mvc_mid[i]);  // sqrt(a)*cos[theta]
        mvc_mid[i] = mvc[i] * tan(mvc_mid[i]);                              // sqrt(a)*cos[theta]*tan[theta]
    }
    return 0;
}

int samp_mvc_gauss(num_real* mvc, num_real variance, int fdim) {
    plFunction();
    rand_gaussian(mvc, 2 * fdim, sqrt(variance));
    return 0;
}

int samp_mvc_sphere(num_real* mvc, num_real Rc2, int fdim) {
    plFunction();
    rand_sphere(mvc, 2 * fdim, sqrt(Rc2));
    return 0;
}


int solve_transform(                                                 ///< solve transform problem
    num_complex* H, num_complex* dH, num_complex* ddH,               // backups
    num_complex* S, num_real* L, num_complex* dL, num_complex* ddL,  // eigenproblem in representation::adiabatic
    num_real* T, num_real* E, num_real* dE, num_real* ddE,           // eigenproblem in representation::diabatic
    num_real* V, num_real* dV, num_real* ddV,                        // original representation::diabatic
    num_real* nr, num_real* np, num_real* nm,                        // phase space
    int rep_type, int level,                                         // flags
    int rdim, int fdim,                                              // shapes
    num_real* workr, num_complex* workc,                             // temperoray array
    bool refered) {
    plFunction();
    int FF = fdim * fdim, NFF = rdim * FF, NNFF = rdim * NFF;
    bool ifbackH = false;

    switch (rep_type) {
        case representation::diabatic: {
            if (fdim < 0) {
                num_real mixangle = -0.5f * atan(2.0f * V[1] / (V[0] - V[3]));
                if (mixangle < 0) mixangle = mixangle + phys::math::halfpi;
                T[0] = cos(mixangle);
                T[1] = sin(mixangle);
                T[2] = -T[1];
                T[3] = T[0];
                E[0] = T[0] * T[0] * V[0] + T[1] * T[1] * V[3] - 2.0f * T[0] * T[1] * V[1];
                E[1] = T[1] * T[1] * V[0] + T[0] * T[0] * V[3] + 2.0f * T[0] * T[1] * V[1];
            } else {
                EigenSolve(E, T, V, fdim);  // sign/sort of eigen problem is not so important
            }
            if (ifbackH) {  ///< @warning record in H, dH, ddH, not recommended!!
                for (int i = 0; i < FF; ++i) H[i] = phys::math::iu * V[i];
                if (level > 0)
                    for (int i = 0; i < NFF; ++i) dH[i] = phys::math::iu * dV[i];
                if (level > 1)
                    for (int i = 0; i < NNFF; ++i) ddH[i] = phys::math::iu * ddV[i];
            }
            break;
        }

        case representation::adiabatic: {
            if (fdim < 0) {  // @warning
                num_real mixangle = -0.5f * atan(2.0f * V[1] / (V[0] - V[3]));
                if (mixangle < 0) mixangle = mixangle + phys::math::halfpi;
                T[0] = cos(mixangle);
                T[1] = sin(mixangle);
                T[2] = -T[1];
                T[3] = T[0];
                E[0] = T[0] * T[0] * V[0] + T[1] * T[1] * V[3] - 2.0f * T[0] * T[1] * V[1];
                E[1] = T[1] * T[1] * V[0] + T[0] * T[0] * V[3] + 2.0f * T[0] * T[1] * V[1];
            } else {
                double* Told   = workr;
                double* TtTold = workr + FF;

                for (int i = 0; i < FF; ++i) Told[i] = T[i];  // backup old T matrix
                EigenSolve(E, T, V, fdim);                    // solve new eigen problem
                if (refered) {                                ///< refer the sign and order of the previous step

                    // calculate permutation matrix = rountint(T^ * Told)
                    ARRAY_MATMUL_TRANS1(TtTold, T, Told, fdim, fdim, fdim);

                    double vset = 0.1 * std::sqrt(1.0e0 / fdim);
                    for (int i = 0; i < fdim; ++i) {
                        double maxnorm = 0;
                        int csr1 = 0, csr2 = 0, csr12 = 0;
                        for (int k1 = 0, k1k2 = 0; k1 < fdim; ++k1) {
                            for (int k2 = 0; k2 < fdim; ++k2, ++k1k2) {
                                if (std::abs(TtTold[k1k2]) > maxnorm) {  // vmax must be larger than sqrt(1/fdim)
                                    maxnorm = std::abs(TtTold[k1k2]);
                                    csr1 = k1, csr2 = k2, csr12 = k1k2;
                                }
                            }
                        }
                        double vsign = copysign(1.0f, TtTold[csr12]);
                        for (int k2 = 0, k1k2 = csr1 * fdim; k2 < fdim; ++k2, ++k1k2) TtTold[k1k2] = 0;
                        for (int k1 = 0, k1k2 = csr2; k1 < fdim; ++k1, k1k2 += fdim) TtTold[k1k2] = 0;
                        TtTold[csr12] = vsign * vset;
                    }
                    for (int i = 0; i < FF; ++i) TtTold[i] = round(TtTold[i] / vset);

                    // // calculate permutation matrix = rountint(T^ * Told)
                    // ARRAY_MATMUL_TRANS1(TtTold, T, Told, fdim, fdim, fdim);
                    // for (int i = 0; i < fdim; ++i) {  // for corresponding to diabatic
                    //     int kcsr   = 0;
                    //     double tmp = std::abs(TtTold[i * fdim]);
                    //     for (int k = 0; k < fdim; ++k) {
                    //         if (std::abs(TtTold[i * fdim + k]) > tmp) {
                    //             tmp  = std::abs(TtTold[i * fdim + k]);
                    //             kcsr = k;
                    //         }
                    //     }
                    //     for (int k = 0; k < fdim; ++k) {
                    //         TtTold[i * fdim + k] = (k == kcsr) ? copysign(1.0f, TtTold[i * fdim + k]) : 0.0f;
                    //     }
                    // }

                    // adjust order of eigenvectors
                    double* Tnew = workr;
                    ARRAY_MATMUL(Tnew, T, TtTold, fdim, fdim, fdim);
                    for (int i = 0; i < FF; ++i) T[i] = Tnew[i];

                    // adjust order of eigenvalues
                    for (int i = 0; i < FF; ++i) TtTold[i] = std::abs(TtTold[i]);
                    double* Enew = workr;
                    ARRAY_MATMUL(Enew, E, TtTold, 1, fdim, fdim);
                    for (int i = 0; i < fdim; ++i) E[i] = Enew[i];
                }
            }

            // calc H = E - im*nacv*np / nm in first step
            for (int i = 0, idx = 0; i < fdim; ++i)
                for (int j = 0; j < fdim; ++j, ++idx) H[idx] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);

            num_real *dEi = dE, *dVi = dV;
            for (int i = 0; i < rdim; ++i) {
                // calc dEi = T^*dVi*T
                num_real* Matr_tmp = workr;  // used as temperoray array
                ARRAY_MATMUL(Matr_tmp, dVi, T, fdim, fdim, fdim);
                ARRAY_MATMUL_TRANS1(dEi, T, Matr_tmp, fdim, fdim, fdim);

                // calc H = E - im*nacv*np/nm in second step
                num_real vi = np[i] / nm[i];
                for (int j = 0; j < fdim; ++j) {  // @TODO optimize
                    for (int k = 0; k < fdim; ++k) {
                        if (j == k) continue;
                        H[j * fdim + k] -= phys::math::im * vi * dEi[j * fdim + k] / (E[k] - E[j]);
                    }
                }
                dVi += FF, dEi += FF;
            }
            EigenSolve(L, S, H, fdim);

            // if necessary, calc ddE
            if (level > 1)
                for (int i = 0; i < rdim; ++i) {
                    num_real *ddVij = ddV + i * NFF, *ddVji = ddV + i * FF;
                    num_real *ddEij = ddE + i * NFF, *ddEji = ddE + i * FF;
                    for (int j = i; j < rdim; ++i) {
                        num_real* Matr_tmp = workr;  // used as temperoray array
                        ARRAY_MATMUL(Matr_tmp, ddVij, T, fdim, fdim, fdim);
                        ARRAY_MATMUL_TRANS1(ddEij, T, Matr_tmp, fdim, fdim, fdim);
                        if (i != j)
                            for (int k = 0; k < FF; ++k) ddEji[k] = ddEij[k];
                        ddVij += FF;
                        ddVji += NFF;
                        ddEij += FF;
                        ddEji += NFF;
                    }
                }
            if (ifbackH) {  // record in H, dH, ddH
                if (level > 0)
                    for (int i = 0; i < NFF; ++i) dH[i] = phys::math::iu * dE[i];
                if (level > 1)
                    for (int i = 0; i < NNFF; ++i) ddH[i] = phys::math::iu * ddE[i];
            }
            break;
        }
        case representation::onthefly: {
            // @TODO track sign of nacv

            // calc H = E - im*nacv*np/nm in first step
            for (int i = 0, idx = 0; i < fdim; ++i)
                for (int j = 0; j < fdim; ++j, ++idx) H[idx] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);
            num_real* dEi = dE;
            for (int i = 0; i < rdim; ++i) {
                // calc H = E - im*nacv*np/nm in second step
                num_real vi = np[i] / nm[i];
                for (int j = 0; j < fdim; ++j) {  // little slow, it doesn't a matter
                    for (int k = 0; k < fdim; ++k) {
                        if (j == k) continue;
                        H[j * fdim + k] -= phys::math::im * vi * dEi[j * fdim + k] / (E[k] - E[j]);
                    }
                }
                dEi += FF;
            }
            EigenSolve(L, S, H, fdim);
            if (level > 1) LOG(FATAL) << "not support!";
            if (ifbackH) {  // record in H, dH, ddH
                if (level > 0)
                    for (int i = 0; i < NFF; ++i) dH[i] = phys::math::iu * dE[i];
                if (level > 1)
                    for (int i = 0; i < NNFF; ++i) ddH[i] = phys::math::iu * ddE[i];
            }
            break;
        }
        case representation::force:
            break;
        case representation::density:
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}


/*
    [ref]: J Phys Chem Lett. 2020 Oct 1;11(19):8283-8291. doi: 10.1021/acs.jpclett.0c02533.
    this transform solves eigen solution for (phase corrected) effective Hamiltonian.
*/
int solve_transform_correctphase(                                    ///< solve transform problem with pc
    num_complex* H, num_complex* rho,                                // additional variables
    num_complex* S, num_real* L, num_complex* dL, num_complex* ddL,  // eigenproblem in representation::adiabatic
    num_real* T, num_real* E, num_real* dE, num_real* ddE,           // eigenproblem in representation::diabatic
    num_real* V, num_real* dV, num_real* ddV,                        // original representation::diabatic
    num_real* nr, num_real* np, num_real* nm,                        // phase space
    int rep_type, int level,                                         // flags
    int rdim, int fdim,                                              // shapes
    num_real* workr, num_complex* workc                              // temperoray array
) {
    plFunction();

    int FF = fdim * fdim, NFF = rdim * FF, NNFF = rdim * NFF;
    switch (rep_type) {
        case representation::adiabatic: {
            for (int i = 0; i < FF; ++i) workr[i] = T[i];
            EigenSolve(E, T, V, fdim);  // @TODO SORT ORDER @Problem
            for (int i = 0; i < fdim; ++i) {
                double ovlp = 0.0f;
                for (int j = 0; j < fdim; ++j) ovlp += workr[j * fdim + i] * T[j * fdim + i];
                if (ovlp < 0.0f)
                    for (int j = 0; j < fdim; ++j) T[j * fdim + i] = -T[j * fdim + i];
            }

            num_complex* Heff = workc;  // as temperoray array

            // calc H = E - im*nacv*np/nm in first step
            for (int i = 0, idx = 0; i < fdim; ++i)
                for (int j = 0; j < fdim; ++j, ++idx) H[idx] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);
            num_real *dEi = dE, *dVi = dV;
            for (int i = 0; i < rdim; ++i) {
                // calc dEi = T^*dVi*T
                num_real* Matr_tmp = workr;  // used as temperoray array
                ARRAY_MATMUL(Matr_tmp, dVi, T, fdim, fdim, fdim);
                ARRAY_MATMUL_TRANS1(dEi, T, Matr_tmp, fdim, fdim, fdim);
                num_real vi = np[i] / nm[i];
                for (int j = 0; j < fdim; ++j) {  // little slow, it doesn't a matter
                    for (int k = 0; k < fdim; ++k) {
                        if (j == k) continue;
                        H[j * fdim + k] -= phys::math::im * vi * dEi[j * fdim + k] / (E[k] - E[j]);
                    }
                }
                dVi += FF, dEi += FF;
            }
            // calc phase corrected Hamiltonian
            double Ekin = 0.0f;
            for (int j = 0; j < rdim; ++j) Ekin += 0.5f * np[j] * np[j] / nm[j];
            double Eavg = REAL_OF(ARRAY_TRACE2(rho, H, fdim, fdim));
            for (int i = 0, idxH = 0; i < fdim; ++i) {
                // if not possible just set all zero
                double s2 = 1.0 + (Eavg - E[i]) / Ekin;
                double s  = (s2 > 0) ? sqrt(s2) : 10e-12;
                for (int k = 0; k < fdim; ++k, ++idxH) {
                    Heff[idxH] = (i == k) ? (H[idxH] - (E[i] + 2 * Ekin * s)) : H[idxH];
                }
            }
            EigenSolve(L, S, Heff, fdim);
            break;
        }
        case representation::onthefly: {
            num_complex* Heff = workc;  // as temperoray array

            // calc H = E - im*nacv*np/nm in first step
            for (int i = 0, idx = 0; i < fdim; ++i)
                for (int j = 0; j < fdim; ++j, ++idx) H[idx] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);
            num_real* dEi = dE;
            for (int i = 0; i < rdim; ++i) {
                // calc H = E - im*nacv*np/nm in second step
                num_real vi = np[i] / nm[i];
                for (int j = 0; j < fdim; ++j) {  // little slow, it doesn't a matter
                    for (int k = 0; k < fdim; ++k) {
                        if (j == k) continue;
                        H[j * fdim + k] -= phys::math::im * vi * dEi[j * fdim + k] / (E[k] - E[j]);
                    }
                }
                dEi += FF;
            }
            // calc phase corrected Hamiltonian
            double Ekin = 0.0f;
            for (int j = 0; j < rdim; ++j) Ekin += 0.5f * np[j] * np[j] / nm[j];
            double Eavg = REAL_OF(ARRAY_TRACE2(rho, H, fdim, fdim));  // @TODO:CHECK?
            for (int i = 0, idxH = 0; i < fdim; ++i) {
                // if not possible just set all zero
                double s2 = 1.0 + (Eavg - E[i]) / Ekin;
                double s  = (s2 > 0) ? sqrt(s2) : 10e-12;
                for (int k = 0; k < fdim; ++k, ++idxH) {
                    Heff[idxH] = (i == k) ? (H[idxH] - (E[i] + 2 * Ekin * s)) : H[idxH];
                }
            }
            EigenSolve(L, S, Heff, fdim);
            break;
        }
        default:
            LOG(FATAL);
    }
    return 0;
}

int solve_Ut(num_complex* U,                      ///< time-evolve propagator
             num_complex* S, num_real* L,         // eigen solution in representation::adiabatic
             num_real* T, num_real* E,            // eigen solution in representation::diabatic
             num_real dtime,                      // time step
             int rep_type,                        // represation
             int rdim, int fdim,                  // shapes
             num_real* workr, num_complex* workc  // temperoray array
) {
    plFunction();

    int FF = fdim * fdim, NFF = rdim * FF, NNFF = rdim * NFF;
    switch (rep_type) {
        case representation::diabatic: {
            // calculate U
            num_complex* invexpiEdt = workc;
            for (int i = 0; i < fdim; ++i)
                invexpiEdt[i] = cos(E[i] * dtime) - phys::math::im * sin(E[i] * dtime);  // invexpiEdt = exp(-im*E*dt)
            ARRAY_MATMUL3_TRANS2(U, T, invexpiEdt, T, fdim, fdim, 0, fdim);
            break;
        }
        case representation::adiabatic:
        case representation::onthefly: {
            // calculate U
            num_complex* invexpiLdt = workc;  // @WARNNING
            for (int i = 0; i < fdim; ++i)
                invexpiLdt[i] = cos(L[i] * dtime) - phys::math::im * sin(L[i] * dtime);  // invexpiLdt = exp(-im*L*dt)
            ARRAY_MATMUL3_TRANS2(U, S, invexpiLdt, S, fdim, fdim, 0, fdim);              // @TODO check
            break;
        }
        default:  // representation::force, representation::density
            LOG(FATAL);
    }
    return 0;
}

int update_eac(num_complex* eac, num_complex* U, int rdim, int fdim, num_real* workr, num_complex* workc) {
    plFunction();

    ARRAY_MATMUL(workc, U, eac, fdim, fdim, 1);
    for (int i = 0; i < fdim; ++i) eac[i] = workc[i];
    return 0;
}

// !> evolution of gamma matrix during dt
int update_rho(num_complex* rho, num_complex* U, int rdim, int fdim, num_real* workr, num_complex* workc) {
    plFunction();

    ARRAY_MATMUL(workc, U, rho, fdim, fdim, fdim);
    ARRAY_MATMUL_TRANS2(rho, workc, U, fdim, fdim, fdim);
    return 0;
}

// !> difference in evolution
int update_drho(num_complex* drho, num_complex* rho, num_complex* U, int rdim, int fdim, num_real* workr,
                num_complex* workc) {
    plFunction();

    ARRAY_MATMUL(workc, U, rho, fdim, fdim, fdim);
    ARRAY_MATMUL_TRANS2(drho, workc, U, fdim, fdim, fdim);
    for (int i = 0; i < fdim * fdim; ++i) drho[i] -= rho[i];
    return 0;
}

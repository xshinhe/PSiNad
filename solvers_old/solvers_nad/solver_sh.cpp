#include "solver_sh.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

namespace hopping {

int solve_Ut(num_complex* U,               ///< time-evolve propagator
             num_complex* S, num_real* L,  // eigen solution in representation::adiabatic
             const num_real& dtime,        // time step
             const int& fdim,
             num_complex* workc  // temperoray array
) {
    int FF = fdim * fdim;
    // calculate U
    num_complex* invexpiLdt = workc;  // @WARNNING
    for (int i = 0; i < fdim; ++i) {
        invexpiLdt[i] = cos(L[i] * dtime) - phys::math::im * sin(L[i] * dtime);  // invexpiLdt = exp(-im*L*dt)
    }
    ARRAY_MATMUL3_TRANS2(U, S, invexpiLdt, S, fdim, fdim, 0, fdim);
    return 0;
}

};  // namespace hopping

int EigenSolve_Refered(num_real* E, num_real* T, num_real* V, num_real* workr, const int& fdim) {
    plFunction();

    int FF           = fdim * fdim;
    num_real* Told   = workr;
    num_real* TtTold = workr + FF;

    for (int i = 0; i < FF; ++i) Told[i] = T[i];  // backup old T matrix
    EigenSolve(E, T, V, fdim);                    // solve eigen problem

    ARRAY_MATMUL_TRANS1(TtTold, T, Told, fdim, fdim, fdim);  // permutation matrix = rount(T^ * Told)

    // updated algorithm
    num_real vset = 0.1 * std::sqrt(1.0e0 / fdim);
    for (int i = 0; i < fdim; ++i) {
        num_real maxnorm = 0;
        int csr1 = 0, csr2 = 0, csr12 = 0;
        for (int k1 = 0, k1k2 = 0; k1 < fdim; ++k1) {
            for (int k2 = 0; k2 < fdim; ++k2, ++k1k2) {
                if (std::abs(TtTold[k1k2]) > maxnorm) {  // vmax must be larger than sqrt(1/fdim)
                    maxnorm = std::abs(TtTold[k1k2]);
                    csr1 = k1, csr2 = k2, csr12 = k1k2;
                }
            }
        }
        num_real vsign = copysign(1.0f, TtTold[csr12]);
        for (int k2 = 0, k1k2 = csr1 * fdim; k2 < fdim; ++k2, ++k1k2) TtTold[k1k2] = 0;
        for (int k1 = 0, k1k2 = csr2; k1 < fdim; ++k1, k1k2 += fdim) TtTold[k1k2] = 0;
        TtTold[csr12] = vsign * vset;
    }
    for (int i = 0; i < FF; ++i) TtTold[i] = round(TtTold[i] / vset);

    /*
    for (int i = 0; i < fdim; ++i) {  // loop in i-col vector
        // find max absvalue of i-col vector
        int kcsr   = 0;
        num_real tmp = std::abs(TtTold[i * fdim + 0]);
        for (int k = 0; k < fdim; ++k) {
            if (std::abs(TtTold[i * fdim + k]) > tmp) {
                tmp  = std::abs(TtTold[i * fdim + k]);
                kcsr = k;
            }
        }
        // round i-col vector to {0,...,0,1,....}
        for (int k = 0; k < fdim; ++k) {
            TtTold[i * fdim + k] = (k == kcsr) ? copysign(1.0f, TtTold[i * fdim + k]) : 0.0f;
        }
    }
    */

    // exchange the order of eigenvectors
    num_real* Tnew = workr;
    ARRAY_MATMUL(Tnew, T, TtTold, fdim, fdim, fdim);
    for (int i = 0; i < FF; ++i) T[i] = Tnew[i];

    // exchange the order of eigenvalues
    for (int i = 0; i < FF; ++i) TtTold[i] = std::abs(TtTold[i]);
    num_real* Enew = workr;
    ARRAY_MATMUL(Enew, E, TtTold, 1, fdim, fdim);
    for (int i = 0; i < fdim; ++i) E[i] = Enew[i];

    return 0;
}


Hopping_Solver::Hopping_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    // default FSSH surface hopping
    Fadd1 = F + 1;

    std::string type = Param_GetT(std::string, parm, "type", "fssh");
    hopping_strategy = hopping::_dict.at(type);

    Param_GetV(reflect, parm, false);
    Param_GetV(pcorrect, parm, false);
    Param_GetV(terminate, parm, false);
    Param_GetV(dish_prefer_deh1, parm, false);
    Param_GetV(dish_prefer_deh2, parm, true);
    Param_GetV(dishw, parm, 1.0f);
    Param_Reset(rep_type, representation::adiabatic);

    eom_type = elec_eom::rho;

    switch (hopping_strategy) {
        case hopping::FSSH:
            LOG(WARNING) << "FSSH";
            break;
        case hopping::GFSH:
            LOG(WARNING) << "GFSH";
            LOG(FATAL);
            break;
        case hopping::DISH:
            LOG(WARNING) << "DISH";
            break;
        case hopping::NDM:
            LOG(WARNING) << "NDM";
            LOG(FATAL);
            break;
        case hopping::SCDM:
            LOG(WARNING) << "SCDM";
            LOG(FATAL);
            break;
        case hopping::AFSSH:
            LOG(WARNING) << "AFSSH";
            break;
    }

    ALLOCATE_PTR_TO_VECTOR(time_coh, FF);
    ALLOCATE_PTR_TO_VECTOR(time_los, FF);
    ALLOCATE_PTR_TO_VECTOR(arg_nr, NFF);
    ALLOCATE_PTR_TO_VECTOR(arg_np, NFF);
    ALLOCATE_PTR_TO_VECTOR(direction, N);
    ALLOCATE_PTR_TO_VECTOR(drho_diss, FF);
    ALLOCATE_PTR_TO_VECTOR(Uh, FF);
    ALLOCATE_PTR_TO_VECTOR(dnp_diss, N);
    ALLOCATE_PTR_TO_VECTOR(probs_dish, F);
    ALLOCATE_PTR_TO_VECTOR(idx_dish, F);

    std::string suffix = type;
    save               = save + "_" + name() + suffix + "_" + pForceField->tag;
};

Hopping_Solver::~Hopping_Solver(){};

int Hopping_Solver::init_occ2eac(const int& itraj) {
    for (int i = 0; i < F; ++i) mvc[i] = (i == occ0) ? 1 : 0;  // first sampling
    samp_mvc_focus(mvc, F);
    eac_mvc(eac0, mvc, F);
    return 0;
};

int Hopping_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);

    memset(arg_nr, 0, NFF * sizeof(num_complex));
    memset(arg_np, 0, NFF * sizeof(num_complex));
    memset(time_coh, 0, FF * sizeof(num_real));
    memset(time_los, 0, FF * sizeof(num_real));
    memset(probs_dish, 0, F * sizeof(num_real));
    memset(idx_dish, 0, F * sizeof(int));

    rho_eac(rho, eac, F);
    if (ini_type == elec_init::d2a) {
        num_real rand_tmp;
        num_real sum = 0.0f;
        rand_uniform(&rand_tmp);
        for (int i = 0; i < F; ++i) {
            sum += NORM_OF(eac[i]);
            if (rand_tmp < sum) {
                occt = i;  // project from diabatic occ0 to adiabatic occt
                break;
            }
        }
    } else {
        occt = occ0;
    }
    return 0;
};

int Hopping_Solver::SE_Hamiltonian(const bool& refered) {  //@done
    plFunction();

    // @todo: if refered?
    if (rep_type == representation::adiabatic) {
        if (refered) {
            EigenSolve_Refered(E, T, V, workr, F);
        } else {
            EigenSolve(E, T, V, F);
        }
    }

    // calc H = E - im*nacv*np/nm in first step
    for (int i = 0, ij = 0; i < F; ++i) {
        for (int j = 0; j < F; ++j, ++ij) {
            H[ij] = ((i == j) ? phys::math::iu * E[i] : phys::math::iz);  // copy eigenvalues
        }
    }

    if (rep_type == representation::adiabatic && N > 100 && F < 5) {  // @opted: faster for large N and small F
        // using EigMXd = Eigen::Matrix<num_real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        // Eigen::Map<EigMXd> MapdV(dV, N * F, F);
        // Eigen::Map<EigMXd> MapdE(dE, N * F, F);
        // Eigen::Map<EigMXd> MapdVT1(workr, F, N * F);
        // Eigen::Map<EigMXd> MapdVT2(workr, F * N, F);
        // Eigen::Map<EigMXd> MapT(T, F, F);
        // MapdVT1 = (MapdV * MapT).transpose();
        // MapdVT2 = (MapdVT2 * MapT).eval();
        // MapdE   = MapdVT1.transpose();

        /// < O(NFFFF) but fast
        num_real* TtxT = workr;
        for (int i1 = 0, i1F = 0, i1k1i2k2 = 0; i1 < F; ++i1, i1F += F) {  // highly optimized
            for (int k1 = 0, k1F = 0; k1 < F; ++k1, k1F += F) {
                for (int i2 = 0, i1i2 = i1F; i2 < F; ++i2, ++i1i2) {
                    for (int k2 = 0, k1k2 = k1F; k2 < F; ++k2, ++k1k2, ++i1k1i2k2) {
                        TtxT[i1k1i2k2] = T[i1i2] * T[k1k2];  // CONJ_OF(T[i1i2]) * T[k1k2]
                    }
                }
            }
        }
        ARRAY_MATMUL(dE, dV, TtxT, N, FF, FF);

        num_real* nvdE = workr;
        for (int i = 0; i < FF; ++i) nvdE[i] = 0;
        for (int j = 0, jk = 0; j < N; ++j) {
            num_real vj = np[j] / nm[j];
            for (int k = 0; k < FF; ++k, ++jk) nvdE[k] += vj * dE[jk];
        }
        for (int i = 0, ik = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++ik) {
                if (i == k) continue;
                H[ik] -= phys::math::im * nvdE[ik] / (E[k] - E[i]);
            }
        }
    } else {  // @general. ~O(NF^3)
        num_real *dEi = dE, *dVi = dV;
        for (int i = 0; i < N; ++i) {
            // for adiabatic rep, calc dE via: dEi = T^*dVi*T
            // but for onthefly rep, dE is previously obtained, so needn't do transform
            if (rep_type == representation::adiabatic) {
                num_real* Matr_tmp = workr;  // used as temperoray array
                ARRAY_MATMUL(Matr_tmp, dVi, T, F, F, F);
                ARRAY_MATMUL_TRANS1(dEi, T, Matr_tmp, F, F, F);
            }

            num_real vi = np[i] / nm[i];
            for (int j = 0, jk = 0; j < F; ++j) {
                for (int k = 0; k < F; ++k, ++jk) {
                    if (j == k) continue;
                    H[jk] -= phys::math::im * vi * dEi[jk] / (E[k] - E[j]);
                }
            }
            dVi += FF, dEi += FF;
        }
    }
    return 0;
}

int Hopping_Solver::phase_correction() {  // @done

    num_real Ekin = nucl_Ekin();
    num_real Epes = REAL_OF(H[occt * Fadd1]);  // phase-correction for surface hopping
    /**
     * @note for phase-correction Ehrenfest:
     *      Epes = REAL_OF(ARRAY_TRACE2(rho, H, F, F));
     */

    for (int i = 0, ii = 0; i < F; ++i, ii += Fadd1) {
        num_real s2 = 1.0 + (Epes - E[i]) / Ekin;
        num_real sc = (s2 > 0) ? sqrt(s2) : 10e-12;  // scaling factor
        H[ii]       = -2 * Ekin * sc;
    }
    return 0;
}

// @brief: generate hopping state from iocc (but with change current state)
int Hopping_Solver::hopping_choose(const int& iocc, num_complex* rhox, num_complex* H, num_real& dt) {
    plFunction();
    int ichoose = -1;
    num_real rand_tmp;
    num_real rhoii = REAL_OF(rhox[iocc * (F + 1)]);
    rand_uniform(&rand_tmp);

    switch (hopping_strategy) {
        // choose with fewest switch algorithm
        case hopping::FSSH:
        case hopping::AFSSH:
        case hopping::NDM:   // @note: not as original NDM, but an improved correction: 10.1063/1.1648306, see at
                             // II.F.3
        case hopping::SCDM:  // SCDM use only SE probability to choose state
        {
            num_real sumprob = 0.0f;
            for (int n = 0; n < F; ++n) {
                num_real prob = (n == iocc) ? 0.0f : -2.0f * IMAG_OF(rhox[n * F + iocc] * H[iocc * F + n]) / rhoii * dt;
                // hopping safely cut-off
                prob = (prob > 1.0f) ? 1.0f : ((prob < 0.0f) ? 0.0f : prob);
                sumprob += prob;
                if (rand_tmp < sumprob) {
                    ichoose = n;
                    break;
                }
            }
            break;
        }
        case hopping::DISH: {
            {
                plScope("dish-time");
                time_calc();  // calculate characteristic time of DISH
            }

            int randi;

            int neff = 0;
            for (int n = 0; n < F; ++n) {
                {
                    plScope("dish-rand");
                    rand_poisson(&randi, 1, time_los[n * Fadd1]);
                }
                // LOG(WARNING) << "poi = " << randi << " ? time " << time_coh[n * Fadd1];
                if (randi <= time_coh[n * Fadd1]) {
                    probs_dish[n]    = REAL_OF(rho[n * Fadd1]);
                    idx_dish[neff++] = n;
                }
            }

            if (neff > 0) {  // choose more than 1

                // LOG(WARNING);

                rand_catalog(&randi, 1, true, 0, neff - 1);
                int decoh_state = idx_dish[randi];

                if (rand_tmp < probs_dish[decoh_state]) {  // project in decoh_state
                    if (decoh_state == occt) {             // needn't hop, directly collapse
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik) {
                                rho[ik] = (i == decoh_state && k == decoh_state) ? phys::math::iu : phys::math::iz;
                            }
                        }
                        ichoose = -1;
                    } else {  // otherwise, try hop and then collapse
                        ichoose           = decoh_state;
                        dish_project_self = false;
                    }
                } else {  // it will project out decoh_state

                    num_real norm_pop = (1 - REAL_OF(rho[decoh_state * Fadd1]));

                    if (decoh_state != occt) {  // projected out decoh_state, and needn't hop
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik) {
                                rho[ik] = (i == decoh_state || k == decoh_state) ? phys::math::iz : rho[ik] / norm_pop;
                            }
                        }
                        ichoose = -1;
                    } else {  // project out occt, need hop to k

                        rand_uniform(&rand_tmp);
                        rand_tmp *= norm_pop;

                        num_real sumprob = 0.0f;
                        for (int k = 0; k < F; ++k) {
                            if (k == decoh_state) continue;
                            sumprob += REAL_OF(rho[k * Fadd1]);
                            if (rand_tmp <= sumprob) {
                                ichoose           = k;  // project out
                                dish_project_self = true;
                                break;
                            }
                        }
                    }
                }

                time_coh[decoh_state * Fadd1] = 0.0f;

            }  // otherwise ichoose = -1
            for (int i = 0; i < FF; ++i) time_coh[i] += dt;
            break;
        }
        default: {
            LOG(FATAL);
        }
    }
    return ichoose;
}

int Hopping_Solver::hopping_direction(num_real* direction, const int& to) {
    plFunction();

    int from = occt;
    switch (hopping_strategy) {
        case hopping::FSSH:
        case hopping::DISH: {
            for (int i = 0; i < N; ++i) direction[i] = dE[i * FF + from * F + to];
            break;
        }
        case hopping::AFSSH: {
            for (int i = 0; i < N; ++i) {  //
                direction[i] = REAL_OF(arg_np[i * FF + to * Fadd1] - arg_np[i * FF + from * Fadd1]);
            }
            break;
        }
        case hopping::NDM: {
            for (int i = 0; i < N; ++i) direction[i] = np[i];
            break;
        }
        case hopping::SCDM: {
            // 1) direction set as unit vector along to d
            num_real dnorm = 0.0f;
            for (int j = 0, jft = from * F + to; j < N; ++j, jft += FF) {
                direction[j] = dE[jft] / (E[to] - E[from]);
                dnorm += direction[j] * direction[j];
            }
            dnorm = sqrt(dnorm);
            for (int j = 0; j < N; ++j) direction[j] /= dnorm;

            num_real pprojd = 0.0f;
            for (int j = 0; j < N; ++j) {
                pprojd += np[j] * direction[j];  // here direction is unit vector
            }

            // d*a0*P^d * \hat{d} \pm P\hat{P}
            if (pprojd > 0) {
                for (int j = 0; j < N; ++j) direction[j] = dnorm * 1 * pprojd * direction[j] + np[j];
            } else {
                for (int j = 0; j < N; ++j) direction[j] = dnorm * 1 * pprojd * direction[j] - np[j];
            }
            break;
        }
        default: {
            LOG(FATAL);
        }
    }
    return 0;
}

int Hopping_Solver::hopping_impulse(num_real* direction, const int& occ_to) {
    plFunction();

    if (occ_to == occt) return 0;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    num_real coeffa = 0.0f, coeffb = 0.0f, coeffc = E[occ_to] - E[occt];
    for (int i = 0; i < N; ++i) {
        coeffa += 0.5f * direction[i] * direction[i] / nm[i];
        coeffb += np[i] / nm[i] * direction[i];
    }
    // normalizing for safe
    coeffb /= coeffa, coeffc /= coeffa;

    num_real coeffd = coeffb * coeffb - 4 * coeffc;
    if (coeffd > 0) {  // else frustrated hopping, do nothing
        num_real x1 = 0.5f * (-coeffb + sqrt(coeffd)), x2 = 0.5f * (-coeffb - sqrt(coeffd));
        num_real x = (std::abs(x1) < std::abs(x2)) ? x1 : x2;
        for (int i = 0; i < N; ++i) np[i] += x * direction[i];

        occt = occ_to;

        if (hopping_strategy == hopping::AFSSH) {  // AFSSH also impulse for argumented variables
            num_complex *arg_nrj = arg_nr, *arg_npj = arg_np;
            for (int j = 0; j < N; ++j) {
                num_complex subr = arg_nrj[occt * Fadd1];
                num_complex subp = arg_npj[occt * Fadd1];
                for (int i = 0, ii = 0; i < F; ++i, ii += Fadd1) {  //
                    arg_nrj[ii] -= subr, arg_npj[ii] -= subp;
                }
                arg_nrj += FF, arg_npj += FF;
            }
        }

    } else if (reflect) {  // 2008Algorithm
        LOG(FATAL);        // DON NOT USE IT
        num_real x = -coeffb;
        for (int i = 0; i < N; ++i) { np[i] += x * direction[i]; }
    } else {
        // LOG(WARNING) << "reject";
    }  // 1990Algorithm, do nothing
    return 0;
}

int Hopping_Solver::time_calc() {
    plFunction();

    num_real Ekin = nucl_Ekin();
    num_real Epes = E[occt];
    num_real Etot = Ekin + Epes;

    switch (hopping_strategy) {
        case hopping::NDM:
        case hopping::SCDM: {
            // tau_{iK}
            for (int i = 0; i < F; ++i) {
                if (i == occt) continue;
                if (hopping_strategy == hopping::NDM)
                    time_los[i * Fadd1] = 1.0f / std::abs(E[occt] - E[i]) * Etot / Ekin;
                if (hopping_strategy == hopping::SCDM)
                    time_los[i * Fadd1] = 1.0f / std::abs(E[occt] - E[i]) + 4.0f / Ekin;
            }
            // occt
            time_los[occt] = 0.0f;
            for (int i = 0; i < F; ++i) {
                if (i == occt) continue;
                time_los[occt * Fadd1] -= REAL_OF(rho[i * Fadd1] / rho[occt * Fadd1]) / time_los[i * Fadd1];
            }
            time_los[occt * Fadd1] = 1.0f / time_los[occt * Fadd1];
            // (i,k) i!=k
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) {
                    if (i == k) continue;
                    time_los[ik] =
                        2 * time_los[i * Fadd1] * time_los[k * Fadd1] / (time_los[i * Fadd1] + time_los[k * Fadd1]);
                }
            }
            break;
        }
        case hopping::DISH: {
            for (int i = 0, ik = 0; i < F; ++i) {
                time_los[i * Fadd1] = 0.0f;
                for (int k = 0; k < F; ++k, ++ik) {
                    if (i == k) continue;

                    // using tau_{ik} as: doi:10.1063/1.470177, you must init mod_W[] at first
                    time_los[ik] = 0.0f;
                    for (int j = 0; j < N; ++j) {
                        time_los[ik] += std::abs(dE[j * FF + i * Fadd1] - dE[j * FF + k * Fadd1]) /
                                        std::abs(np[j] * dishw * dishw / phys::math::twopi);
                    }
                    time_los[ik] = 1.0f / time_los[ik];

                    time_los[i * Fadd1] += REAL_OF(rho[k * Fadd1]) * 1.0f / time_los[ik];
                }
                time_los[i * Fadd1] = 1.0f / time_los[i * Fadd1];
            }
            // for safe
            for (int i = 0; i < FF; ++i) {
                time_los[i] = std::abs(time_los[i]);
                if (time_los[i] > 1.0e8) time_los[i] = 1.0e8;
            }
            break;
        }
        default: {
            LOG(FATAL);
        }
    }
    return 0;
}

int Hopping_Solver::coherent_evolve(const num_real& dt) {
    plFunction();

    // assuming hamiltonian is already built & solved
    // update rho
    hopping::solve_Ut(U, S, L, dt, F, workc);
    update_rho(rho, U, N, F, workr, workc);
    // update_eac(eac, U, N, F, workr, workc);

    switch (hopping_strategy) {
        case hopping::AFSSH: {  // unique to AFSSH

            hopping::solve_Ut(Uh, S, L, 0.5e0 * dt, F, workc);

            // using velocity-verlet to update argumented variables
            num_real* dEj        = dE;
            num_complex *arg_nrj = arg_nr, *arg_npj = arg_np;
            for (int j = 0; j < N; ++j) {
                update_rho(arg_nrj, Uh, N, F, workr, workc);
                for (int i = 0; i < FF; ++i) arg_nrj[i] += arg_npj[i] / nm[j] * 0.5e0 * dt;

                num_real dEjoo = dEj[occt * Fadd1];

                update_rho(arg_npj, U, N, F, workr, workc);
                num_complex* FM1 = workc;
                num_complex* FM2 = workc + FF;
                ARRAY_MATMUL(FM1, rho, dEj, F, F, F);
                ARRAY_MATMUL(FM2, dEj, rho, F, F, F);
                for (int i = 0; i < FF; ++i) arg_npj[i] -= (0.5e0 * (FM1[i] + FM2[i]) - dEjoo * rho[i]) * dt;

                update_rho(arg_nrj, Uh, N, F, workr, workc);
                for (int i = 0; i < FF; ++i) arg_nrj[i] += arg_npj[i] / nm[j] * 0.5e0 * dt;

                dEj += FF, arg_nrj += FF, arg_npj += FF;
            }
        }
        default: {
        }
    }
    return 0;
}

int Hopping_Solver::decoherent_evolve(const num_real& dt) {
    plFunction();

    switch (hopping_strategy) {
        case hopping::NDM:
        case hopping::SCDM: {  // unique to NDM & SCDM

            plScope("xxx");

            time_calc();
            for (int i = 0; i < FF; ++i) drho_diss[i] = -rho[i] / time_los[i] * dt;

            memset(dnp_diss, 0, N * sizeof(num_real));
            for (int k = 0; k < F; ++k) {
                if (k == occt) continue;

                hopping_direction(direction, k);
                num_real vdotsk = 0.0f;
                for (int j = 0; j < N; ++j) vdotsk += np[j] / nm[j] * direction[j];

                num_real sum = REAL_OF(rho[k * Fadd1]) / time_los[occt * F + k] * E[occt] +
                               REAL_OF(drho_diss[k * Fadd1]) / dt * E[k];

                // no-contribution, since for adiabatic representation E
                // sum += REAL_OF(drho_diss[occt * F + k]) * H[occt * F + k];
                // for (int kp = 0; kp < F; ++kp) {
                //     if (kp == k) continue;
                //     sum += REAL_OF(drho_diss[k * F + kp]) * H[k * F + kp];
                // }

                for (int j = 0; j < N; ++j) dnp_diss[j] = -sum * direction[j] / vdotsk;
            }
            for (int i = 0; i < FF; ++i) rho[i] += drho_diss[i];
            for (int i = 0; i < N; ++i) np[i] += dnp_diss[i];
            break;
        }
        default: {
            break;
        }
    }
    return 0;
}

int Hopping_Solver::hopping_collapse(const num_real& dt, const int& to) {
    plFunction();

    switch (hopping_strategy) {
        case hopping::AFSSH: {  // hopping & collapsing are seperated (randomly generated)

            num_real rand_tmp, sumprob = 0.0f;
            rand_uniform(&rand_tmp);
            num_real* dEj        = dE;
            num_complex *arg_nrj = arg_nr, *arg_npj = arg_np;
            for (int j = 0; j < N; ++j) {
                for (int n = 0; n < F; ++n) {
                    num_real sign = copysign(1.0f, REAL_OF(arg_nrj[0] - arg_nrj[3]) * REAL_OF(arg_npj[0] - arg_npj[3]));
                    num_real prob = 0.5e0 * sign * (-dEj[0] + dEj[3]) * REAL_OF(arg_nrj[0] - arg_nrj[3]) * dt;
                    // if (prob > 0) ??
                    sumprob += prob;
                }
                dEj += FF, arg_nrj += FF, arg_npj += FF;
            }
            if (rand_tmp < sumprob) {  // do collapse
                for (int i = 0, ik = 0; i < F; ++i) {
                    for (int k = 0; k < F; ++k, ++ik) {
                        rho[ik] = (i == occt && k == occt) ? phys::math::iu : phys::math::iz;
                    }
                }
                memset(arg_nr, 0, NFF * sizeof(num_complex));
                memset(arg_np, 0, NFF * sizeof(num_complex));
            }
            break;
        }
        case hopping::DISH: {  // collapsing is same with hopping
            {
                if (!dish_project_self) {
                    if (dish_prefer_deh1 ||
                        occt == to) {  // if you prefer it, or succeed in nucl hopping event, then do elec project
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik)
                                rho[ik] = (i == to && k == to) ? phys::math::iu : phys::math::iz;
                        }
                    } else {
                        num_real norm_pop = 1.0f - REAL_OF(rho[to * Fadd1]);
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik)
                                rho[ik] = (i == to || k == to) ? phys::math::iz : rho[ik] / norm_pop;
                        }
                    }
                } else {
                    if (dish_prefer_deh2 || occt == to) {  // then keep self in occt
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik)
                                rho[ik] = (i == to && k == to) ? phys::math::iu : phys::math::iz;
                        }
                    } else {
                        for (int i = 0, ik = 0; i < F; ++i) {
                            for (int k = 0; k < F; ++k, ++ik)
                                rho[ik] = (i == occt && k == occt) ? phys::math::iu : phys::math::iz;
                        }
                    }
                }
            }
            break;
        }
        default: {  // no operations
        }
    }
    return 0;
}

int Hopping_Solver::ff_calc1(const int& level,
                             const bool& refered) {  // hamiltonian calculated for electronic DOFs
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

    SE_Hamiltonian(refered);
    if (pcorrect) phase_correction();
    return 0;
}

int Hopping_Solver::kernel_fssh(num_complex* rhox, const int& F) {
    switch (hopping_strategy) {
        case hopping::DISH:
            for (int i = 0; i < FF; ++i) rhox[i] = rho[i];
            break;
        default: {
            for (int i = 0, ik = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++ik) {
                    rhox[ik] = (i == k) ? ((i == occt) ? phys::math::iu : phys::math::iz) : rho[ik];
                }
            }
        }
    }
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    return 0;
};

int Hopping_Solver::check_break(int& succ) {
    plFunction();
    if (succ != 0 || !ARRAY_ISFINITE(nr, N) || !ARRAY_ISFINITE(np, N) || !ARRAY_ISFINITE(rho, F * F)) { succ = -1; }
    return succ;
}

int Hopping_Solver::kernel0(num_complex* rhox, const int& F) {
    for (int i = 0, idx = 0; i < F; ++i) {  // @todo: bug
        for (int j = 0; j < F; ++j, ++idx) { rhox[idx] = ((i == j) ? phys::math::iu : phys::math::iz); }
    }
    return 0;
};

int Hopping_Solver::kernelt(num_complex* rhox, const int& F) { return kernel_fssh(rhox, F); };

int Hopping_Solver::ff_calc2() {
    plFunction();

    for (int j = 0, joo = occt * (F + 1); j < N; ++j, joo += FF) fmean[j] = dE[joo];
    for (int j = 0; j < N; ++j) nf[j] = grad[j] + fmean[j];
    return 0;
};

int Hopping_Solver::traj(NAD_TCFer& tcfer, const int& N, const int& F) {
    init(itraj);
    ff_calc1(level), ff_calc2();

    int hop_count   = 0;
    num_real thres  = 0.75f * std::abs(V[0] - V[3]);
    num_real delta0 = 1e10, delta1 = 1e10, delta2 = 1e10;

    int succ = sampler(0, tcfer);

    for (int istep_dummy = 0; istep_dummy < nstep; ++istep_dummy) {
        if (succ == 0) istep = istep_dummy;  // otherwise froze istep

        bool if_small   = (delta1 < thres);
        int nsec        = if_small ? nrespa : 1;
        num_real dtime  = if_small ? smalldt : dt;
        num_real dtimeh = if_small ? 0.5f * smalldt : 0.5f * dt;

        for (int isec = 0; isec < nsec; ++isec) {
            if (succ == 0) succ = update_p(dtimeh);
            if (succ == 0) succ = update_r(dtimeh);
            if (succ == 0 && dyn_type == -1 && pThermo->dothermo(istep_dummy)) succ = update_thermo(dtime);
            if (succ == 0) succ = update_r(dtimeh);
            if (succ == 0) succ = ff_calc1(level, true);
            delta0 = delta1, delta1 = delta2, delta2 = std::abs(V[0] - V[3]);

            if (succ == 0) EigenSolve(L, S, H, F);
            if (succ == 0) coherent_evolve(dtime);
            if (succ == 0) decoherent_evolve(dtime);

            if (succ == 0) {
                int to = hopping_choose(occt, rho, H, dtime);
                if (to >= 0) {  // meaningful choose
                    hopping_direction(direction, to);
                    hopping_impulse(direction, to);
                    hopping_collapse(dtime, to);
                }
            }
            if (succ == 0) succ = ff_calc2();
            if (succ == 0) succ = update_p(dtimeh);
        }

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

        // for time distribution
        if (terminate && std::abs(nr[0]) > 10.0f) {
            if (succ == 0) succ = sampler(1, tcfer);
            if (!ofs_SAMP.is_open()) LOG(FATAL) << "ofs_SAMP should be open";

            ofs_SAMP << FMT(0) << hop_count << FMT(0) << itraj;
            ofs_SAMP << FMT(12) << (istep + 1) * dt  //
                     << FMT(12) << nr[0]             //
                     << FMT(12) << np[0]             //
                     << FMT(12) << nr0[0]            //
                     << FMT(12) << np0[0];
            for (int i = 0; i < tcfer.lentcf; ++i)
                ofs_SAMP << FMT(12) << REAL_OF(tcfer.val[i]) << FMT(12) << IMAG_OF(tcfer.val[i]);
            ofs_SAMP << std::endl;
            break;
        }
    }

    final(itraj);
    return 0;
}
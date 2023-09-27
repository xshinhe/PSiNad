#include "mespimd_solver.h"

#include "../utils/definitions.h"

/**
 * @brief implement of Multi-Electronic-State Path-Integral-Molecular-Dynamics (MES-PIMD)
 * @reference
 *  1)  J. Chem. Phys. 148, 102319 (2018). doi:10.1063/1.5005059.
 *  2)  Chin. J. Chem. Phys. 31, 4, 446-456 (2018). doi: 10.1063/1674-0068/31/cjcp1805122
 */

using namespace ARRAY_EG;

MESPIMDTraj_Solver::MESPIMDTraj_Solver(const Param& iparm, Model* pM) : PIMDTraj_Solver(iparm, pM) {
    pForceField = dynamic_cast<Nad_ForceField*>(pM);

    factor_flag = Param_GetT(std::string, parm, "factor_flag", "first");
    factor_type = MESPIMDFactorizePolicy::_dict.at(factor_flag);

    save = name() + "_" + pForceField->tag;

    // get dimensions
    F    = pForceField->get_F();
    FF   = F * F;
    NF   = N * F;
    NFF  = N * F * F;
    NNF  = N * N * F;
    NNFF = N * N * F * F;
    FFFF = F * F * F * F;

    ALLOCATE_PTR_TO_VECTOR(Vs, P * FF);
    ALLOCATE_PTR_TO_VECTOR(dVs, P * NFF);
    // ALLOCATE_PTR_TO_VECTOR(ddVs, P * NNFF);

    ALLOCATE_PTR_TO_VECTOR(Es, P * F);
    ALLOCATE_PTR_TO_VECTOR(Ts, P * FF);
    ALLOCATE_PTR_TO_VECTOR(dEs, P * NFF);

    ALLOCATE_PTR_TO_VECTOR(O, P * FF);
    ALLOCATE_PTR_TO_VECTOR(OO, P * FF);
    ALLOCATE_PTR_TO_VECTOR(OOb, P * FF);
    ALLOCATE_PTR_TO_VECTOR(OOe, P * FF);
    ALLOCATE_PTR_TO_VECTOR(OObe, PP * FF);
    ALLOCATE_PTR_TO_VECTOR(dO, PN * FF);
    ALLOCATE_PTR_TO_VECTOR(V2, P * FF);
    ALLOCATE_PTR_TO_VECTOR(hRdOO, P * FF);
    ALLOCATE_PTR_TO_VECTOR(hRcdOO, P * FF);
    ALLOCATE_PTR_TO_VECTOR(hRdOVO, P * FF);
    ALLOCATE_PTR_TO_VECTOR(hRcdOVO, P * FF);

    ALLOCATE_PTR_TO_VECTOR(O_pos, P * FF);
    ALLOCATE_PTR_TO_VECTOR(dV_pos, PN * FF);

    ALLOCATE_PTR_TO_VECTOR(rho, FF);
    ALLOCATE_PTR_TO_VECTOR(rho_op, FF);
    ALLOCATE_PTR_TO_VECTOR(mat_fD, FF);
    ALLOCATE_PTR_TO_VECTOR(mat_fA, FF);
    ALLOCATE_PTR_TO_VECTOR(eac0, F);
    ALLOCATE_PTR_TO_VECTOR(rho0, FF);
}

MESPIMDTraj_Solver::~MESPIMDTraj_Solver(){};

int MESPIMDTraj_Solver::init(const int& itraj) {
    PIMDTraj_Solver::init(itraj);  // invalid sampling from BO_ForceField's ForceField_init

    for (int i = 0, ibead_j = 0; i < P; ++i) {  ///< Nad_ForceField's ForceField_init
        pForceField->ForceField_init(nr0, np0, nm, rho0, eac0, occ0, N, F, itraj);
        for (int j = 0; j < N; ++j, ++ibead_j) {
            nrs[ibead_j] = nr0[j], nps[ibead_j] = np0[j], nms[ibead_j] = masswgt[i] * nm[j];
        }
    }

    if (FLAGS_r != "") rst_read(itraj);
    return 0;
}

/**
 * @brief solve transform on i-th bead
 *
 * @param Oi_pos positive-part extract from \prod exp(-0.5*betap*V)
 * @param Oi approximate to exp(-0.5*betap*Vi)
 * @param dOi approximate to d_R exp(-0.5*betap*Vi)
 * @param Vi potential matrix (FF) on i-th bead
 * @param dVi force tensor (NFF) on i-th bead
 * @return success flag
 */
int MESPIMDTraj_Solver::ff_Oi(int i, num_real* Oi_pos, num_real* dVi_pos, num_real* Oi, num_real* dOi, num_real* Vi,
                              num_real* dVi) {
    switch (factor_type) {
        case MESPIMDFactorizePolicy::Diag: {
            /**
             * diagonal expansion to O = exp(-0.5 * betap * V)
             *  O_pos = exp(-0.5 * betap * Ediag)
             *  O     = T* exp(-0.5 * betap * Ediag)*T^
             *  dO    = dT* exp(-0.5 * betap * Ediag)*T^
             *          + T* dexp(-0.5 * betap * Ediag)*T^
             *          + T* exp(-0.5 * betap * Ediag)*dT^
             *          = Ttau* exp(-0.5 * betap * Ediag)*T^
             *          + T* dexp(-0.5 * betap * Ediag)*T^
             *          + T* exp(-0.5 * betap * Ediag)*tau^T^
             */
            Eigen::Map<EigMXr> Map_Oi_pos(Oi_pos, F, F);
            Eigen::Map<EigMXr> Map_Oi(Oi, F, F);
            Eigen::Map<EigMXr> Map_Ei(Es + i * F, F, 1);
            Eigen::Map<EigMXr> Map_Ti(Ts + i * FF, F, F);

            num_real maxE = (Map_Ei.array()).abs().maxCoeff();

            EigMXr MatE = Map_Ei.asDiagonal();
            EigMXr expd = (-0.5f * betap * Map_Ei.array()).exp().matrix();
            EigMXr Expd = expd.asDiagonal();

            Map_Oi_pos = Expd;
            Map_Oi     = Map_Ti * Expd * Map_Ti.transpose();

            int L = 10, M = 20;
            for (int j = 0, jFF = 0; j < N; ++j, jFF += FF) {
                Eigen::Map<EigMXr> Map_dVij_pos(dVi_pos + jFF, F, F);
                Eigen::Map<EigMXr> Map_dVij(dVi + jFF, F, F);
                Eigen::Map<EigMXr> Map_dOij(dOi + jFF, F, F);
                auto Mat_dEij = Map_Ti.transpose() * Map_dVij * Map_Ti;

                num_real betap_small = betap;
                for (int i = 0; i < L; ++i) betap_small /= 2;

                Map_dOij           = EigMXr::Zero(F, F);
                num_real prefactor = 1;
                for (int m = 1; m < M; ++m) {
                    EigMXr Mtmp = EigMXr::Zero(F, F);

                    prefactor *= (-0.5e0 * betap_small) / m;
                    if (abs(prefactor) * maxE < 1.0e-14) break;

                    for (j = 0; j < m; ++j) {
                        EigMXr Epow1 = (Map_Ei.array().pow(j)).matrix().asDiagonal();
                        EigMXr Epow2 = (Map_Ei.array().pow(m - j - 1)).matrix().asDiagonal();
                        Mtmp += Epow1 * Mat_dEij * Epow2;
                    }
                    Map_dOij += prefactor * Mtmp;
                }

                for (int i = 0; i < L; ++i) {
                    EigMXr expd_small = (-0.5f * betap_small * Map_Ei.array()).exp().matrix();
                    EigMXr Expd_small = expd_small.asDiagonal();
                    Map_dOij          = (Map_dOij * Expd_small + Expd_small * Map_dOij).eval();
                    betap_small *= 2;
                }
                Map_dOij     = (Map_Ti * Map_dOij * Map_Ti.transpose()).eval();
                Map_dVij_pos = (Mat_dEij.diagonal()).asDiagonal();
            }
            break;
        }
        case MESPIMDFactorizePolicy::Hyper: {
            /**
             * first-order expansion to O = exp(-0.5 * betap * V)
             *  O_pos = exp(-0.5 * betap * Vdiag)
             *  O     = (\prod_{j<i} W^{ij}) * exp(-0.5 * betap * Vdiag)
             *  dO    = -0.5 * betap * [dVoffd + (1 - 0.5 * betap * Voffd) * dVdiag]
             *          * exp(-0.5 * betap * Vdiag)
             */

            Eigen::Map<EigMXr> Map_Oi_pos(Oi_pos, F, F);
            Eigen::Map<EigMXr> Map_Oi(Oi, F, F);
            Eigen::Map<EigMXr> Map_Vi(Vi, F, F);

            EigMXr diag = Map_Vi.diagonal();  // diagonal vector
            EigMXr Diag = diag.asDiagonal();  // diagonal matrix
            auto Offd   = Map_Vi - Diag;      // off-diagonal matrix

            EigMXr expd = (-0.5f * betap * diag.array()).exp().matrix();
            EigMXr Expd = expd.asDiagonal();

            EigMXr Expo  = EigMXr::Identity(F, F);
            EigMXr W22a  = EigMXr::Identity(2, 2);
            EigMXr W22b  = EigMXr::Identity(2, 2);
            EigMXr Mcosh = (((-0.5e0 * betap * Map_Vi).array()).cosh()).matrix();
            EigMXr Msinh = (((-0.5e0 * betap * Map_Vi).array()).sinh()).matrix();
            for (int i = 0; i < F; ++i) {
                for (int k = 0; k < i; ++k) {
                    W22a(0, 0) = Expo(i, i);
                    W22a(1, 1) = Expo(k, k);
                    W22a(0, 1) = Expo(i, k);
                    W22a(1, 0) = Expo(k, i);
                    W22b(0, 0) = W22b(1, 1) = Mcosh(i, k);
                    W22b(0, 1) = W22b(1, 0) = Msinh(i, k);
                    W22a                    = (W22a * W22b).eval();
                    Expo(i, i)              = W22a(0, 0);
                    Expo(k, k)              = W22a(1, 1);
                    Expo(i, k)              = W22a(0, 1);
                    Expo(k, i)              = W22a(1, 0);
                }
            }

            Map_Oi_pos = Expd;
            Map_Oi     = Expo * Expd;

            for (int j = 0, jFF = 0; j < N; ++j, jFF += FF) {
                Eigen::Map<EigMXr> Map_dVij_pos(dVi_pos + jFF, F, F);
                Eigen::Map<EigMXr> Map_dVij(dVi + jFF, F, F);
                Eigen::Map<EigMXr> Map_dOij(dOi + jFF, F, F);

                EigMXr diag_dVij = Map_dVij.diagonal();     // diagonal vector
                EigMXr Diag_dVij = diag_dVij.asDiagonal();  // diagonal matrix

                auto dexpd   = (expd.array() * diag_dVij.array()).matrix();
                EigMXr dExpd = dexpd.asDiagonal();

                EigMXr dExpo = EigMXr::Zero(F, F);
                for (int i = 0; i < F; ++i) {
                    for (int k = 0; k < i; ++k) {
                        EigMXr dikExpo = EigMXr::Identity(F, F);
                        for (int i2 = 0; i2 < F; ++i2) {
                            for (int k2 = 0; k2 < i2; ++k2) {
                                W22a(0, 0) = dikExpo(i2, i2);
                                W22a(1, 1) = dikExpo(k2, k2);
                                W22a(0, 1) = dikExpo(i2, k2);
                                W22a(1, 0) = dikExpo(k2, i2);
                                if (i == i2 && k == k2) {
                                    W22b(0, 0) = W22b(1, 1) = Mcosh(i2, k2) * Map_dVij(i2, k2);
                                    W22b(0, 1) = W22b(1, 0) = Msinh(i2, k2) * Map_dVij(i2, k2);
                                } else {
                                    W22b(0, 0) = W22b(1, 1) = Msinh(i2, k2);
                                    W22b(0, 1) = W22b(1, 0) = Mcosh(i2, k2);
                                }
                                W22a            = (W22a * W22b).eval();
                                dikExpo(i2, i2) = W22a(0, 0);
                                dikExpo(k2, k2) = W22a(1, 1);
                                dikExpo(i2, k2) = W22a(0, 1);
                                dikExpo(k2, i2) = W22a(1, 0);
                            }
                        }
                        dExpo += dikExpo;
                    }
                }
                Map_dVij_pos = Diag_dVij;
                Map_dOij     = -0.5f * betap * (dExpo * Expd + Expo * dExpd);
            }
            break;
        }
        default:
        case MESPIMDFactorizePolicy::First: {
            /**
             * first-order expansion to O = exp(-0.5 * betap * V)
             *  O_pos = exp(-0.5 * betap * Vdiag)
             *  O     = (1 - 0.5 * betap * Voffd) * exp(-0.5 * betap * Vdiag)
             *  dO    = -0.5 * betap * [dVoffd + (1 - 0.5 * betap * Voffd) * dVdiag]
             *          * exp(-0.5 * betap * Vdiag)
             */

            Eigen::Map<EigMXr> Map_Oi_pos(Oi_pos, F, F);
            Eigen::Map<EigMXr> Map_Oi(Oi, F, F);
            Eigen::Map<EigMXr> Map_Vi(Vi, F, F);

            EigMXr diag = Map_Vi.diagonal();  // diagonal vector
            EigMXr Diag = diag.asDiagonal();  // diagonal matrix
            auto Offd   = Map_Vi - Diag;      // off-diagonal matrix

            EigMXr expd = (-0.5f * betap * diag.array()).exp().matrix();
            EigMXr Expd = expd.asDiagonal();
            auto Expo   = EigMXr::Identity(F, F) - 0.5f * betap * Offd;

            Map_Oi_pos = Expd;
            Map_Oi     = Expo * Expd;

            for (int j = 0, jFF = 0; j < N; ++j, jFF += FF) {
                Eigen::Map<EigMXr> Map_dVij_pos(dVi_pos + jFF, F, F);
                Eigen::Map<EigMXr> Map_dVij(dVi + jFF, F, F);
                Eigen::Map<EigMXr> Map_dOij(dOi + jFF, F, F);

                EigMXr diag_dVij = Map_dVij.diagonal();     // diagonal vector
                EigMXr Diag_dVij = diag_dVij.asDiagonal();  // diagonal matrix
                auto Offd_dVij   = Map_dVij - Diag_dVij;    // off-diagonal matrix

                auto dexpd   = (expd.array() * diag_dVij.array()).matrix();
                EigMXr dExpd = dexpd.asDiagonal();

                Map_dVij_pos = Diag_dVij;
                Map_dOij     = -0.5f * betap * (Map_dVij * Expd + Expo * dExpd);
            }
            break;
        }
    }
    return 0;
}

int MESPIMDTraj_Solver::ff_OO() {
    Eigen::Map<EigMXr> Map_fD(mat_fD, F, F);
    Eigen::Map<EigMXr> Map_fA(mat_fA, F, F);
    Map_fD         = EigMXr::Identity(F, F);
    Map_fA         = EigMXr::Identity(F, F);
    EigMXr OOprods = EigMXr::Identity(F, F);
    for (int i = 0, J = 0, idxV = 0, idxdV = 0; i < P; ++i, idxV += FF, idxdV += NFF) {
        num_real* Oi_pos  = O_pos + idxV;
        num_real* Oi      = O + idxV;
        num_real* OOi     = OO + idxV;
        num_real* OObi    = OOb + idxV;
        num_real* Vi      = Vs + idxV;
        num_real* V2i     = V2 + idxV;
        num_real* dOi     = dO + idxdV;
        num_real* dVi_pos = dV_pos + idxdV;
        num_real* dVi     = dVs + idxdV;
        Eigen::Map<EigMXr> Map_Oi_pos(Oi_pos, F, F);
        Eigen::Map<EigMXr> Map_Oi(Oi, F, F);
        Eigen::Map<EigMXr> Map_OOi(OOi, F, F);
        Eigen::Map<EigMXr> Map_OObi(OObi, F, F);
        Eigen::Map<EigMXr> Map_Vi(Vi, F, F);
        Eigen::Map<EigMXr> Map_V2i(V2i, F, F);

        ff_Oi(i, Oi_pos, dVi_pos, Oi, dOi, Vi, dVi);

        // calculate related quantities
        Map_V2i = Map_Vi * Map_Vi;
        Map_OOi = Map_Oi.transpose() * Map_Oi;
        if (i == 0) {
            OOprods = EigMXr::Identity(F, F);
        } else {
            Eigen::Map<EigMXr> Map_OOlast(OOi - FF, F, F);
            OOprods = (OOprods * Map_OOlast).eval();
        }
        Map_OObi = OOprods;

        ///< calculate: fD & fA
        Map_fD = (Map_fD * Map_Oi_pos.transpose() * Map_Oi_pos).eval();
        Map_fA = (Map_fA * Map_OOi).eval();

        ///< calculate hRdOO & hRdOVO
        Eigen::Map<EigMXr> Map_hRdOOi(hRdOO + idxV, F, F);
        Eigen::Map<EigMXr> Map_hRdOVOi(hRdOVO + idxV, F, F);
        Eigen::Map<EigMXr> Map_hRcdOOi(hRcdOO + idxV, F, F);
        Eigen::Map<EigMXr> Map_hRcdOVOi(hRcdOVO + idxV, F, F);
        Map_hRdOOi   = EigMXr::Zero(F, F);
        Map_hRdOVOi  = EigMXr::Zero(F, F);
        Map_hRcdOOi  = EigMXr::Zero(F, F);
        Map_hRcdOVOi = EigMXr::Zero(F, F);
        for (int j = 0, jFF = 0; j < N; ++j, ++J, jFF += FF) {
            Eigen::Map<EigMXr> Map_dVij(dVi + jFF, F, F);
            Eigen::Map<EigMXr> Map_dOij(dOi + jFF, F, F);

            Map_hRdOOi += 0.5f * (nrs[J]) *                    //
                          (                                    //
                              Map_Oi.transpose() * Map_dOij    //
                              + Map_dOij.transpose() * Map_Oi  //
                          );
            Map_hRdOVOi += 0.5f * (nrs[J]) *                             //
                           (                                             //
                               Map_Oi.transpose() * Map_Vi * Map_dOij    //
                               + Map_dOij.transpose() * Map_Vi * Map_Oi  //
                               + Map_Oi.transpose() * Map_dVij * Map_Oi  //
                           );
            Map_hRcdOOi += 0.5f * (nrs[J] - nr[j]) *            //
                           (                                    //
                               Map_Oi.transpose() * Map_dOij    //
                               + Map_dOij.transpose() * Map_Oi  //
                           );
            Map_hRcdOVOi += 0.5f * (nrs[J] - nr[j]) *                     //
                            (                                             //
                                Map_Oi.transpose() * Map_Vi * Map_dOij    //
                                + Map_dOij.transpose() * Map_Vi * Map_Oi  //
                                + Map_Oi.transpose() * Map_dVij * Map_Oi  //
                            );
        }
    }
    tr_fD     = Map_fD.trace();
    tr_fA     = Map_fA.trace();
    rw_factor = tr_fA / tr_fD;

    /* for OOe */
    for (int i = P - 1, iFF = i * FF; i >= 0; --i, iFF -= FF) {
        Eigen::Map<EigMXr> Map_OOei(OOe + iFF, F, F);
        if (i == P - 1) {
            OOprods = EigMXr::Identity(F, F);
        } else {
            Eigen::Map<EigMXr> Map_OOlast(OO + iFF + FF, F, F);
            OOprods = (Map_OOlast * OOprods).eval();
        }
        Map_OOei = OOprods;
    }

    /* for OObe */
    const int PFF     = P * FF;
    const int Padd1FF = PFF + FF;
    for (int d = 0, dFF = 0; d < P; ++d, dFF += FF) {  // @highly optimized
        for (int i = 0, iFF = 0, iPadd1FF = 0; i < P - d; ++i, iFF += FF, iPadd1FF += Padd1FF) {
            Eigen::Map<EigMXr> Map_OObeid(OObe + iPadd1FF + dFF, F, F);  ///< only use half of saves
            if (d <= 1) {
                Map_OObeid = EigMXr::Identity(F, F);
            } else if (d == 2) {
                Eigen::Map<EigMXr> Map_OO1(OO + iFF + FF, F, F);
                Map_OObeid = Map_OO1;
            } else {
                Eigen::Map<EigMXr> Map_OObe1(OObe + iPadd1FF + FF + FF, F, F);
                Eigen::Map<EigMXr> Map_OObe2(OObe + iPadd1FF + PFF + dFF, F, F);
                Map_OObeid = Map_OObe1 * Map_OObe2;
            }
        }
    }
    return 0;
}

int MESPIMDTraj_Solver::ff_calc1(const int& level) {
    all_K2X();

    // calc centroid r
    ARRAY_CLEAR(nr, N);
    for (int i = 0, ibead_j = 0; i < P; ++i)
        for (int j = 0; j < N; ++j, ++ibead_j) nr[j] += nrs[ibead_j];
    for (int j = 0; j < N; ++j) nr[j] /= P;

    // calc ForceField
    for (int i = 0, idxR = 0, idxV = 0, idxdV = 0; i < P; ++i, idxR += N, idxV += FF, idxdV += NFF) {
        num_real* nri   = nrs + idxR;
        num_real* npi   = nps + idxR;
        num_real* vpesi = vpeses + i;
        num_real* gradi = grads + idxR;
        num_real* hessi = hesses + 0;  ///< @todo don't use hessian
        num_real* Vi    = Vs + idxV;
        num_real* dVi   = dVs + idxdV;
        // ddV: a nullptr

        pForceField->ForceField_npes(vpesi, gradi, hessi, nri, npi, level, N, itraj, istep);
        pForceField->ForceField_epes(Vi, dVi, ddVs, nri, level, N, F, itraj, istep);

        if (factor_type == MESPIMDFactorizePolicy::Diag) {
            num_real* Ei  = Es + i * F;
            num_real* Ti  = Ts + i * FF;
            num_real* dEi = dEs + i * NFF;
            EigenSolve(Ei, Ti, Vi, F);
        }
    }

    ff_OO();

    Eigen::Map<EigMXr> Map_fD(mat_fD, F, F);
    for (int i = 0, idxR = 0, idxdVj = 0, J = 0; i < P; ++i, idxR += N) {
        num_real* gradi = grads + idxR;
        num_real* nfi   = nfs + idxR;
        for (int j = 0; j < N; ++j, ++J, idxdVj += FF) {
            Eigen::Map<EigMXr> Map_dVij_pos(dV_pos + idxdVj, F, F);
            nfs[J] = (Map_dVij_pos * Map_fD).trace() / (tr_fD);
        }
        for (int j = 0; j < N; ++j) nfi[j] += gradi[j];  ///< add nuclear part and electronic part together
    }
    all_FX2FK();  ///< transfrom f to fk
    return 0;
}

// Q is [P, F, F] dimensional array for unfixed, while [F, F] for fixed flag
num_real MESPIMDTraj_Solver::esti_term1(num_real* Q, const bool& dressed, const bool& fixed) {
    num_real esti = 0.0f;
    for (int i = 0, idxV = 0; i < P; ++i, idxV += FF) {
        Eigen::Map<EigMXr> Map_Q(Q + ((fixed) ? 0 : idxV), F, F);
        Eigen::Map<EigMXr> Map_O(O + idxV, F, F);
        Eigen::Map<EigMXr> Map_OOb(OOb + idxV, F, F);
        Eigen::Map<EigMXr> Map_OOe(OOe + idxV, F, F);
        EigMXr tmpm = Map_Q;  // @bug don't use auto, it will give reference
        if (dressed) tmpm = Map_O.transpose() * Map_Q * Map_O;
        esti += (Map_OOb * tmpm * Map_OOe).trace();
    }
    return esti;
}

// Q1 & Q2 is [P, F, F] dimensional array for unfixed, while [F, F] for fixed flag
num_real MESPIMDTraj_Solver::esti_term2(num_real* Q1, num_real* Q2, const bool& dressed1, const bool& dressed2,
                                        const bool& fixed1, const bool& fixed2) {
    num_real esti = 0.0f;
    for (int i = 0, idxV1 = 0, I = 0; i < P; ++i, idxV1 += FF) {
        Eigen::Map<EigMXr> Map_Q1(Q1 + ((fixed1) ? 0 : idxV1), F, F);
        Eigen::Map<EigMXr> Map_O1(O + idxV1, F, F);
        Eigen::Map<EigMXr> Map_OOb1(OOb + idxV1, F, F);
        Eigen::Map<EigMXr> Map_OOe1(OOe + idxV1, F, F);

        EigMXr tmpm1 = Map_Q1;  // @bug don't use auto, it will give reference
        if (dressed1) tmpm1 = Map_O1.transpose() * Map_Q1 * Map_O1;

        for (int k = 0, idxV2 = 0; k < P; ++k, idxV2 += FF, ++I) {
            Eigen::Map<EigMXr> Map_Q2(Q2 + ((fixed2) ? 0 : idxV2), F, F);
            Eigen::Map<EigMXr> Map_O2(O + idxV2, F, F);
            Eigen::Map<EigMXr> Map_OOb2(OOb + idxV2, F, F);
            Eigen::Map<EigMXr> Map_OOe2(OOe + idxV2, F, F);

            EigMXr tmpm2 = Map_Q2;  // @bug don't use auto, it will give reference
            if (dressed2) tmpm2 = Map_O2.transpose() * Map_Q2 * Map_O2;

            if (i < k) {  // i<k
                Eigen::Map<EigMXr> Map_OObe(OObe + (i * P + k) * FF, F, F);
                esti += (Map_OOb1 * tmpm1 * Map_OObe * tmpm2 * Map_OOe2).trace();
            } else if (i > k) {  // i>k
                Eigen::Map<EigMXr> Map_OObe(OObe + (k * P + i) * FF, F, F);
                esti += (Map_OOb2 * tmpm2 * Map_OObe * tmpm1 * Map_OOe1).trace();
            } else {
                // pass when i == k
            }
        }
    }
    return esti;
}

int MESPIMDTraj_Solver::estimator(const int& isamp, TCFnucl& tcfer) {
    ///< estimator for potential
    // nuclear part: Vn = -\partial_\beta \ln f_n
    num_real esti_Vn = 0;
    for (int i = 0; i < P; ++i) esti_Vn += vpeses[i] / P;

    // electronic part: Ve = -\partial_\beta \ln f_e
    num_real esti_Ve = esti_term1(Vs, true) / (P * tr_fA);

    ///< estimator for primitive kinetic energy
    // primitive K (for both nuclear & electronic dofs)
    num_real esti_Kprim = 0.5f * N * P / pThermo->beta;  // normalization factor in PIMD
    for (int i = 0, J = 0, Jp = N; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++J, ++Jp %= PN) {
            esti_Kprim -= 0.5f * bf2 * nm[j] * (nrs[J] - nrs[Jp]) * (nrs[J] - nrs[Jp]);  // using nrs
            // if (i != 0) esti_Kprim -= 0.5f * bf2 * bfwgt[i] * nms[J] * nks[J] * nks[J]; // using nks
        }
    }
    ///< estimator for virtual kinetic energy
    // const part
    num_real esti_Kvir0 = 0.5f * N / pThermo->beta;
    num_real esti_Kvirn = 0.0f;
    // nuclear part: Kvn = -\beta^{-1} \hat \xi \ln f_n = \hat \xi Vn
    for (int i = 0, ibead_j = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++ibead_j) { esti_Kvirn += 0.5f * nrs[ibead_j] * grads[ibead_j] / P; }
    }
    // nuclear part remove centroid
    num_real esti_Kvirn_c = 0.0f;
    for (int i = 0, ibead_j = 0; i < P; ++i) {
        for (int j = 0; j < N; ++j, ++ibead_j) { esti_Kvirn_c += 0.5f * (nrs[ibead_j] - nr[j]) * grads[ibead_j] / P; }
    }
    // electronic part: Kve = -\beta^{-1} \hat \xi f_e
    num_real esti_Kvire   = -1.0f / (pThermo->beta) * esti_term1(hRdOO, false) / tr_fA;
    num_real esti_Kvire_c = -1.0f / (pThermo->beta) * esti_term1(hRcdOO, false) / tr_fA;

    // f_e^{-1} \partial_\beta^2 f_e
    num_real theta12 = (esti_term1(V2, true) + esti_term2(Vs, Vs, true, true)) / (P * P * tr_fA);

    // -f_e^{-1} \beta^{-1}\hat \xi \partial_\beta f_e
    num_real theta345 = -(esti_term1(hRdOVO, false) + esti_term2(Vs, hRdOO, true, false)) / (pThermo->beta * P * tr_fA);
    num_real theta345_c =
        -(esti_term1(hRcdOVO, false) + esti_term2(Vs, hRcdOO, true, false)) / (pThermo->beta * P * tr_fA);

    Eigen::Map<EigMXr> Map_rho(rho, F, F);
    Eigen::Map<EigMXr> Map_rho_op(rho_op, F, F);
    for (int i = 0; i < F; ++i) {
        for (int j = i; j < F; ++j) {
            Map_rho_op = EigMXr::Zero(F, F);
            Map_rho_op(i, j) += 0.5f;
            Map_rho_op(j, i) += 0.5f;
            Map_rho(i, j) = esti_term1(rho_op, true, true) / (P * tr_fA);
            if (i != j) Map_rho(j, i) = Map_rho(i, j);
        }
    }

    if (isamp == 0) {
        ofs_ENER << FMT(4) << "flag" << FMT(4) << "time"                             //
                 << FMT(4) << "P"                                                    // No. of beads
                 << FMT(4) << "N"                                                    // No. of ndofs
                 << FMT(8) << "beta"                                                 //
                 << FMT(8) << "tr_fD"                                                //
                 << FMT(8) << "rw_factor"                                            //
                 << FMT(8) << "esti_Vn"                                              //
                 << FMT(8) << "esti_Ve"                                              //
                 << FMT(8) << "esti_Kp"                                              //
                 << FMT(8) << "esti_Kvn"                                             //
                 << FMT(8) << "esti_Kvn_c"                                           //
                 << FMT(8) << "esti_Kve"                                             //
                 << FMT(8) << "esti_Kve_c"                                           //
                 << FMT(8) << "theta12"                                              //
                 << FMT(8) << "theta345"                                             //
                 << FMT(8) << "theta345_c";                                          //
        for (int i = 0; i < FF; ++i) ofs_ENER << FMT(8) << utils::concat("rho", i);  // density
        ofs_ENER << std::endl;
    }

    num_real sampunit = dt * sstep * iou.time;
    ofs_ENER << FMT(4) << itraj << FMT(8) << isamp * sampunit   //
             << FMT(4) << P                                     // No. of beads
             << FMT(4) << N                                     // No. of ndofs
             << FMT(8) << pThermo->beta                         //
             << FMT(8) << tr_fD                                 //
             << FMT(8) << rw_factor                             //
             << FMT(8) << esti_Vn                               //
             << FMT(8) << esti_Ve                               //
             << FMT(8) << esti_Kprim                            //
             << FMT(8) << esti_Kvirn                            //
             << FMT(8) << esti_Kvirn_c                          //
             << FMT(8) << esti_Kvire                            //
             << FMT(8) << esti_Kvire_c                          //
             << FMT(8) << theta12                               //
             << FMT(8) << theta345                              //
             << FMT(8) << theta345_c;                           //
    for (int i = 0; i < FF; ++i) ofs_ENER << FMT(8) << rho[i];  // density
    ofs_ENER << std::endl;
    return 0;
}

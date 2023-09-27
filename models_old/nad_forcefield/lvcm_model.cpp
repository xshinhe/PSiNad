#include "lvcm_model.h"

#include "../../utils/definitions.h"
#include "hamiltonian_data.h"


LVCM_ForceField::LVCM_ForceField(const Param& iparm, const int& child) : Nad_ForceField(iparm){};  // for child to call
LVCM_ForceField::LVCM_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    // 1) initial system Hamiltonian
    ALLOCATE_PTR_TO_VECTOR(Hsys, FF);
    ALLOCATE_PTR_TO_VECTOR(ECI, F);
    ALLOCATE_PTR_TO_VECTOR(KCI, N * F);

    try {
        Param_GetV(ffflag, parm, "pyr3");
        fftype = LVCMPolicy::_dict.at(ffflag);
    } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }

    std::string Hsys_flag  = Param_GetT(std::string, parm, "Hsys_flag", "inner");
    std::string H_unit_rev = Param_GetT(std::string, parm, "H_unit", "null");
    LOG(WARNING) << "overwrite unit for reading hamiltonian as <H_unit>" << H_unit_rev << std::endl;
    num_real H_unit = (H_unit_rev == "null") ? iou.ener : phys::unitsys::conv(phys::au::unit, H_unit_rev);
    // revise the H_unit if necessary

    switch (fftype) {
        case LVCMPolicy::PYR3: {
            CHECK_EQ(F, 2);
            CHECK_EQ(N, 3);
            double tmp_unit = phys::au_2_ev;


            CHECK_EQ((F + 1) * (N + 1), sizeof(CI_Pyr3_data) / sizeof(CI_Pyr3_data[0]));

            int icnt = 0;
            lcoeff   = CI_Pyr3_data[icnt++] / tmp_unit;
            for (int i = 0; i < F; ++i) ECI[i] = CI_Pyr3_data[icnt++] / tmp_unit;
            for (int j = 0, idxk = 0; j < N; ++j) {
                mod_W[j] = CI_Pyr3_data[icnt++] / tmp_unit;
                for (int i = 0; i < F; ++i, ++idxk) KCI[idxk] = CI_Pyr3_data[icnt++] / tmp_unit;
                mod_sigmaR[j] = sqrt(0.5f / mod_W[j]);
                mod_sigmaP[j] = sqrt(0.5f * mod_W[j]);
                mod_M[j]      = 1.0f;
                mod_R0[j]     = 0.0f;
                mod_P0[j]     = 0.0f;
            }
            break;
        }
        case LVCMPolicy::PYR4: {
            CHECK_EQ(F, 2);
            CHECK_EQ(N, 4);
            double tmp_unit = phys::au_2_ev;
            CHECK_EQ((F + 1) * (N + 1), sizeof(CI_Pyr4_data) / sizeof(CI_Pyr4_data[0]));

            int icnt = 0;
            lcoeff   = CI_Pyr4_data[icnt++] / tmp_unit;
            for (int i = 0; i < F; ++i) ECI[i] = CI_Pyr4_data[icnt++] / tmp_unit;
            for (int j = 0, idxk = 0; j < N; ++j) {
                mod_W[j] = CI_Pyr4_data[icnt++] / tmp_unit;
                for (int i = 0; i < F; ++i, ++idxk) KCI[idxk] = CI_Pyr4_data[icnt++] / tmp_unit;
                mod_sigmaR[j] = sqrt(0.5f / mod_W[j]);
                mod_sigmaP[j] = sqrt(0.5f * mod_W[j]);
                mod_M[j]      = 1.0f;
                mod_R0[j]     = 0.0f;
                mod_P0[j]     = 0.0f;
            }
            break;
        }
        case LVCMPolicy::PYR24: {
            CHECK_EQ(F, 2);
            CHECK_EQ(N, 24);
            double tmp_unit = phys::au_2_ev;
            CHECK_EQ((F + 1) * (N + 1), sizeof(CI_Pyr24_data) / sizeof(CI_Pyr24_data[0]));

            int icnt = 0;
            lcoeff   = CI_Pyr24_data[icnt++] / tmp_unit;
            for (int i = 0; i < F; ++i) ECI[i] = CI_Pyr24_data[icnt++] / tmp_unit;
            for (int j = 0, idxk = 0; j < N; ++j) {
                mod_W[j] = CI_Pyr24_data[icnt++] / tmp_unit;
                for (int i = 0; i < F; ++i, ++idxk) KCI[idxk] = CI_Pyr24_data[icnt++] / tmp_unit;
                mod_sigmaR[j] = sqrt(0.5f / mod_W[j]);
                mod_sigmaP[j] = sqrt(0.5f * mod_W[j]);
                mod_M[j]      = 1.0f;
                mod_R0[j]     = 0.0f;
                mod_P0[j]     = 0.0f;
            }
            break;
        }
        case LVCMPolicy::CRC2: {
            CHECK_EQ(F, 3);
            CHECK_EQ(N, 2);
            double freq[5] = {0.0129f, 0.0129f, 0.0342f, 0.0561f, 0.0561f};
            double reqb[5] = {0.0f, 14.3514f, -9.9699f, -7.0189f, 0.0f};
            double alpw[5] = {0.4501f, 0.4286f, 0.6204f, 0.4535f, 0.5539f};
            for (int j = 0; j < N; ++j) {
                mod_W[j]      = freq[j] / phys::au_2_ev;
                mod_R0[j]     = reqb[j] / sqrt(mod_W[j]);
                mod_P0[j]     = 0.0f;
                mod_sigmaR[j] = alpw[j] / sqrt(mod_W[j]);
                mod_sigmaP[j] = 0.5f * sqrt(mod_W[j]) / alpw[j];
                mod_M[j]      = 1.0f;
            }
            break;
        }
        case LVCMPolicy::CRC5: {
            CHECK_EQ(F, 3);
            CHECK_EQ(N, 5);
            double freq[5] = {0.0129f, 0.0129f, 0.0342f, 0.0561f, 0.0561f};
            double reqb[5] = {0.0f, 14.3514f, -9.9699f, -7.0189f, 0.0f};
            double alpw[5] = {0.4501f, 0.4286f, 0.6204f, 0.4535f, 0.5539f};
            for (int j = 0; j < N; ++j) {
                mod_W[j]      = freq[j] / phys::au_2_ev;
                mod_R0[j]     = reqb[j] / sqrt(mod_W[j]);
                mod_P0[j]     = 0.0f;
                mod_sigmaR[j] = alpw[j] / sqrt(mod_W[j]);
                mod_sigmaP[j] = 0.5f * sqrt(mod_W[j]) / alpw[j];
                mod_M[j]      = 1.0f;
            }
            break;
        }
        case LVCMPolicy::BEN5: {
            CHECK_EQ(F, 3);
            CHECK_EQ(N, 5);
            double tmp_unit = phys::au_2_ev;
            CHECK_EQ((F + 1) * (N + 1), sizeof(CI_BEN5_data) / sizeof(CI_BEN5_data[0]));

            int icnt = 0;
            lcoeff   = CI_BEN5_data[icnt++] / tmp_unit;  // useless
            for (int i = 0; i < F; ++i) ECI[i] = CI_BEN5_data[icnt++] / tmp_unit;
            for (int j = 0, idxk = 0; j < N; ++j) {
                mod_W[j] = CI_BEN5_data[icnt++] / tmp_unit;
                for (int i = 0; i < F; ++i, ++idxk) KCI[idxk] = CI_BEN5_data[icnt++] / tmp_unit;
                mod_sigmaR[j] = sqrt(0.5f / mod_W[j]);
                mod_sigmaP[j] = sqrt(0.5f * mod_W[j]);
                mod_M[j]      = 1.0f;
                mod_R0[j]     = 0.0f;
                mod_P0[j]     = 0.0f;
            }
            break;
        }
        default:
            LOG(FATAL);
    }

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

int LVCM_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                     int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int LVCM_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                     const int& rdim) {
    V[0] = 0.0f;
    for (int j = 0; j < N; ++j) {
        V[0] += 0.5f * (P[j] * P[j] / mod_M[j] + mod_M[j] * mod_W[j] * mod_W[j] * R[j] * R[j]);
    }

    if (flag < 1) return 0;

    for (int j = 0; j < N; ++j) { dV[j] = mod_M[j] * mod_W[j] * mod_W[j] * R[j]; }

    if (flag < 2) return 0;
    for (int j = 0; j < NN; ++j) ddV[j] = 0;
    int add = N + 1;
    for (int j = 0, idx = 0; j < N; ++j, idx += add) ddV[idx] = mod_M[j] * mod_W[j] * mod_W[j];
    return 0;
}

int LVCM_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                     const int& fdim) {
    switch (fftype) {
        case LVCMPolicy::PYR3:
        case LVCMPolicy::PYR4:
        case LVCMPolicy::PYR24:
            ForceField_epes_PYR(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case LVCMPolicy::CRC2:
            ForceField_epes_CrCO5_2(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case LVCMPolicy::CRC5:
            ForceField_epes_CrCO5_5(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case LVCMPolicy::BEN5:
            ForceField_epes_BEN_5(V, dV, ddV, R, flag, rdim, fdim);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int LVCM_ForceField::ForceField_epes_PYR(double* V, double* dV, double* ddV, double* R, const int& flag,
                                         const int& rdim, const int& fdim) {
    for (int i = 0, idx = 0; i < F; ++i)
        for (int k = 0; k < F; ++k, ++idx)
            V[idx] = (i == k) ? ECI[i] : lcoeff * sqrt(mod_W[0]) * R[0];  // first mode in LC-CI
    for (int j = 0, idxk = 0; j < N; ++j) {
        V[0] += KCI[idxk++] * sqrt(mod_W[j]) * R[j];
        V[3] += KCI[idxk++] * sqrt(mod_W[j]) * R[j];
    }
    if (flag < 1) return 0;

    if (first_call) {
        for (int j = 0, idxk = 0, idxdV = 0; j < N; ++j) {
            dV[idxdV++] = KCI[idxk++];
            dV[idxdV++] = (j == 0) ? lcoeff : 0.0f;
            dV[idxdV++] = (j == 0) ? lcoeff : 0.0f;
            dV[idxdV++] = KCI[idxk++];
        }
        for (int j = 0, idxdV = 0; j < N; ++j)
            for (int i = 0; i < FF; ++i, ++idxdV) dV[idxdV] *= sqrt(mod_W[j]);
    }
    first_call = false;

    if (flag < 2) return 0;

    return 0;
}


int LVCM_ForceField::ForceField_epes_CrCO5_2(double* V, double* dV, double* ddV, double* R, const int& flag,
                                             const int& rdim, const int& fdim) {
    const double EE       = 0.0424f / phys::au_2_ev;
    const double EA       = 0.4344f / phys::au_2_ev;
    const double kappa1   = -0.0169f / phys::au_2_ev;
    const double kappa3   = -0.0272f / phys::au_2_ev;
    const double lambda1a = 0.0328f / phys::au_2_ev;
    const double lambda1b = -0.0978f / phys::au_2_ev;
    const double lambda2a = 0.0095f / phys::au_2_ev;
    const double lambda2b = 0.1014f / phys::au_2_ev;

    double x0 = sqrt(mod_W[0]) * R[0];
    double x1 = sqrt(mod_W[1]) * R[1];

    V[0] = EE - lambda1a * x1;
    V[4] = EE + lambda1a * x1;
    V[8] = EA;
    V[1] = lambda1a * x0;
    V[2] = lambda1b * x1;
    V[5] = lambda1b * x0;
    V[3] = V[1];
    V[6] = V[2];
    V[7] = V[5];

    if (flag < 1) return 0;

    if (first_call) {
        for (int j = 0; j < NFF; ++j) dV[j] = 0.0f;

        dV[0 * FF + 1] = lambda1a * sqrt(mod_W[0]);
        dV[0 * FF + 3] = lambda1a * sqrt(mod_W[0]);
        dV[0 * FF + 5] = lambda1b * sqrt(mod_W[0]);
        dV[0 * FF + 7] = lambda1b * sqrt(mod_W[0]);

        dV[1 * FF + 0] = -lambda1a * sqrt(mod_W[1]);
        dV[1 * FF + 4] = lambda1a * sqrt(mod_W[1]);
        dV[1 * FF + 2] = lambda1b * sqrt(mod_W[1]);
        dV[1 * FF + 6] = lambda1b * sqrt(mod_W[1]);
    }
    first_call = false;

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int LVCM_ForceField::ForceField_epes_CrCO5_5(double* V, double* dV, double* ddV, double* R, const int& flag,
                                             const int& rdim, const int& fdim) {
    const double EE       = 0.0424f / phys::au_2_ev;
    const double EA       = 0.4344f / phys::au_2_ev;
    const double kappa1   = -0.0169f / phys::au_2_ev;
    const double kappa3   = -0.0272f / phys::au_2_ev;
    const double lambda1a = 0.0328f / phys::au_2_ev;
    const double lambda1b = -0.0978f / phys::au_2_ev;
    const double lambda2a = 0.0095f / phys::au_2_ev;
    const double lambda2b = 0.1014f / phys::au_2_ev;

    double x0 = sqrt(mod_W[0]) * R[0];
    double x1 = sqrt(mod_W[1]) * R[1];
    double x2 = sqrt(mod_W[2]) * R[2];
    double x3 = sqrt(mod_W[3]) * R[3];
    double x4 = sqrt(mod_W[4]) * R[4];

    V[0] = EE + kappa1 * x2 - lambda1a * x1 - lambda2a * x3;
    V[4] = EE + kappa1 * x2 + lambda1a * x1 + lambda2a * x3;
    V[8] = EA + kappa3 * x2;
    V[1] = lambda1a * x0 + lambda2a * x4;
    V[2] = lambda1b * x1 + lambda2b * x3;
    V[5] = lambda1b * x0 + lambda2b * x4;
    V[3] = V[1];
    V[6] = V[2];
    V[7] = V[5];

    if (flag < 1) return 0;

    if (first_call) {
        for (int j = 0; j < NFF; ++j) dV[j] = 0.0f;

        dV[0 * FF + 1] = lambda1a * sqrt(mod_W[0]);
        dV[0 * FF + 3] = lambda1a * sqrt(mod_W[0]);
        dV[0 * FF + 5] = lambda1b * sqrt(mod_W[0]);
        dV[0 * FF + 7] = lambda1b * sqrt(mod_W[0]);

        dV[1 * FF + 0] = -lambda1a * sqrt(mod_W[1]);
        dV[1 * FF + 4] = lambda1a * sqrt(mod_W[1]);
        dV[1 * FF + 2] = lambda1b * sqrt(mod_W[1]);
        dV[1 * FF + 6] = lambda1b * sqrt(mod_W[1]);

        dV[2 * FF + 0] = kappa1 * sqrt(mod_W[2]);
        dV[2 * FF + 4] = kappa1 * sqrt(mod_W[2]);
        dV[2 * FF + 8] = kappa3 * sqrt(mod_W[2]);


        dV[3 * FF + 0] = -lambda2a * sqrt(mod_W[4]);
        dV[3 * FF + 4] = lambda2a * sqrt(mod_W[4]);
        dV[3 * FF + 2] = lambda2b * sqrt(mod_W[4]);
        dV[3 * FF + 6] = lambda2b * sqrt(mod_W[4]);

        dV[4 * FF + 1] = lambda2a * sqrt(mod_W[3]);
        dV[4 * FF + 3] = lambda2a * sqrt(mod_W[3]);
        dV[4 * FF + 5] = lambda2b * sqrt(mod_W[3]);
        dV[4 * FF + 7] = lambda2b * sqrt(mod_W[3]);
    }
    first_call = false;

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int LVCM_ForceField::ForceField_epes_BEN_5(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double au_2_ev = 27.21138602f;
    const double lc12    = 0.164f / au_2_ev;
    const double lc23    = 0.154f / au_2_ev;

    for (int i = 0, idx = 0; i < F; ++i)
        for (int k = 0; k < F; ++k, ++idx) {
            if (i == k) {
                V[idx] = ECI[i];
            } else if ((i == 0 && k == 1) || (i == 1 & k == 0)) {
                V[idx] = lc12 * sqrt(mod_W[3]) * R[3];
            } else if ((i == 1 && k == 2) || (i == 2 & k == 1)) {
                V[idx] = lc23 * sqrt(mod_W[4]) * R[4];
            } else {
                V[idx] = 0;
            }
        }
    for (int j = 0, idxk = 0; j < N; ++j) {
        for (int i = 0, Fadd1 = F + 1; i < F; ++i, ++idxk) {  //
            V[i * Fadd1] += KCI[idxk] * sqrt(mod_W[j]) * R[j];
        }
    }
    if (flag < 1) return 0;

    if (first_call) {
        for (int j = 0, idxk = 0, idxdV = 0; j < N; ++j) {
            for (int i = 0; i < F; ++i, ++idxk) {
                for (int k = 0; k < F; ++k, ++idxdV) {
                    if (i == k) {
                        dV[idxdV] = KCI[idxk];
                    } else if ((i == 0 && k == 1) || (i == 1 & k == 0)) {
                        dV[idxdV] = (j == 3) ? lc12 : 0;
                    } else if ((i == 1 && k == 2) || (i == 2 & k == 1)) {
                        dV[idxdV] = (j == 4) ? lc23 : 0;
                    } else {
                        dV[idxdV] = 0;
                    }
                }
            }
        }
        for (int j = 0, idxdV = 0; j < N; ++j)
            for (int i = 0; i < FF; ++i, ++idxdV) dV[idxdV] *= sqrt(mod_W[j]);
    }
    first_call = false;
    // LOG(FATAL);

    if (flag < 2) return 0;

    return 0;
}

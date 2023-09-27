#include "pyrincavity_models.h"

#include "../../utils/definitions.h"
#include "hamiltonian_data.h"


PyrCav_ForceField::PyrCav_ForceField(const Param& iparm, const int& child)
    : Nad_ForceField(iparm){};  // for child to call

PyrCav_ForceField::PyrCav_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    // double phys::au_2_ev = 1; LOG(FATAL);

    try {
        Param_GetV(ffflag, parm, "pc2");
        fftype = PyrCavPolicy::_dict.at(ffflag);
    } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }

    Param_GetV(gcoup, parm, 0.24f);
    gcoup /= phys::au_2_ev;

    Param_GetV(wcav, parm, 0.62f);
    wcav /= phys::au_2_ev;

    Param_GetV(Fmol, parm, 3);  // S0, S1, S2
    Param_GetV(Nmod, parm, 3);  // 3 molecule-modes pyrazine model
    CHECK_EQ(Fmol, 3);
    CHECK_EQ((Nmod - 2) * (Nmod - 3) * (Nmod - 4) * (Nmod - 24), 0);

    // 1) initial system Hamiltonian
    ALLOCATE_PTR_TO_VECTOR(ECI, Fmol);
    ALLOCATE_PTR_TO_VECTOR(KCI, Nmod * (Fmol - 1));  // only on S1&S2
    ALLOCATE_PTR_TO_VECTOR(WCI, Nmod);

    const double* H_data;
    switch (Nmod) {
        case 2:
            H_data = CI_Pyr2_data;
            LOG(WARNING) << "3-modes pyrazine model";
            break;
        case 3:
            H_data = CI_Pyr3_data;
            LOG(WARNING) << "3-modes pyrazine model";
            break;
        case 4:
            H_data = CI_Pyr4_data;
            LOG(WARNING) << "4-modes pyrazine model";
            break;
        case 24:
            H_data = CI_Pyr24_data;
            LOG(WARNING) << "24-modes pyrazine model";
            break;
        default:
            LOG(FATAL);
    }

    int icnt = 0;
    lcoeff   = H_data[icnt++] / phys::au_2_ev;
    ECI[0]   = 0.0f;                                                         // S0
    for (int i = 1; i < Fmol; ++i) ECI[i] = H_data[icnt++] / phys::au_2_ev;  // S1, S2
    for (int j = 0, idxk = 0; j < Nmod; ++j) {
        WCI[j] = H_data[icnt++] / phys::au_2_ev;
        for (int i = 1; i < Fmol; ++i, ++idxk) {  // skip S0
            KCI[idxk] = H_data[icnt++] / phys::au_2_ev;
        }
    }
    // update lcoeff
    double mylcoeff = Param_GetV(mylcoeff, iparm, -1.0f);
    if (mylcoeff > 0) lcoeff = mylcoeff / phys::au_2_ev;

    switch (fftype) {
        case PyrCavPolicy::PC1: {
            LOG(FATAL);
            break;
        }
        case PyrCavPolicy::PC2: {
            Param_GetV(Ncav, parm, 2);
            CHECK_EQ(F, Ncav * Fmol);
            CHECK_EQ(N, Nmod);  // molecule-modes
            for (int j = 0; j < N; ++j) {
                mod_W[j]      = WCI[j];
                mod_sigmaR[j] = sqrt(0.5f / mod_W[j]);
                mod_sigmaP[j] = sqrt(0.5f * mod_W[j]);
                mod_M[j]      = 1.0f;
                mod_R0[j]     = 0.0f;
                mod_P0[j]     = 0.0f;
            }
            break;
        }
        case PyrCavPolicy::PC3: {
            CHECK_EQ(F, Fmol);
            CHECK_EQ(N, Nmod + 1);  // cavity-modes & molecule-modes
            for (int j = 0; j < N; ++j) {
                mod_W[j]      = ((j == 0) ? wcav : WCI[j - 1]);
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

int PyrCav_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                       int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int PyrCav_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
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

int PyrCav_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                       const int& fdim) {
    switch (fftype) {
        case PyrCavPolicy::PC1:  // all treated as quantum
            ForceField_epes_PC1(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case PyrCavPolicy::PC2:
            ForceField_epes_PC2(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case PyrCavPolicy::PC3:  // cav-modes & ci-modes treaded as quantum
            ForceField_epes_PC3(V, dV, ddV, R, flag, rdim, fdim);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int PyrCav_ForceField::ForceField_epes_PC1(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    LOG(FATAL);
    return 0;
}

int PyrCav_ForceField::ForceField_epes_PC2(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    /**
     *
     * R[0]     : molecule-modes (ci-mode)
     * R[1...]  : molecule-modes (bath-modes)
     *
     */
    int Fmolm1 = Fmol - 1;  // Fmolm1 == 2

    if (first_call)
        for (int i = 0; i < FF; ++i) V[i] = 0.0f;

    for (int icav = 0, i = 0, idx = 0; icav < Ncav; ++icav) {
        double Ecav = icav * wcav;
        for (int ilev = 0; ilev < Fmol; ++ilev, ++i) {
            for (int kcav = 0, k = 0; kcav < Ncav; ++kcav) {
                for (int klev = 0; klev < Fmol; ++klev, ++k, ++idx) {
                    if (icav == kcav) {  // in the same cavity modes
                        if (ilev == 0 || klev == 0) {
                            if (ilev == klev) V[idx] = Ecav;
                            continue;  // skip S0
                        }
                        if (ilev == klev) {
                            V[idx] = Ecav + ECI[ilev];  // cavity & atomic energy
                            for (int j = 0, jk = ilev - 1; j < N; ++j, jk += Fmolm1) {
                                V[idx] += KCI[jk] * sqrt(mod_W[j]) * R[j];
                            }
                        } else {
                            V[idx] = lcoeff * sqrt(mod_W[0]) * R[0];  // ci-modes
                        }
                    } else {
                        // continue;
                        if (icav == 0 && kcav == 1 && ilev == 2 && klev == 1) V[idx] = gcoup;
                        if (icav == 0 && kcav == 1 && ilev == 1 && klev == 2) V[idx] = gcoup;
                        if (icav == 1 && kcav == 0 && ilev == 1 && klev == 2) V[idx] = gcoup;
                        if (icav == 1 && kcav == 0 && ilev == 2 && klev == 1) V[idx] = gcoup;
                    }
                }
            }
        }
    }

    // ARRAY_COLOR_SHOW(V, F, F, 1e-6);
    // LOG(FATAL);

    if (flag < 1) return 0;

    if (first_call) {
        for (int i = 0; i < NFF; ++i) dV[i] = 0.0f;  // init to zero

        for (int j = 0, idxdV = 0; j < N; ++j) {
            for (int icav = 0, i = 0, idx = 0; icav < Ncav; ++icav) {
                for (int ilev = 0; ilev < Fmol; ++ilev, ++i) {
                    for (int kcav = 0, k = 0; kcav < Ncav; ++kcav) {
                        for (int klev = 0; klev < Fmol; ++klev, ++k, ++idx, ++idxdV) {
                            if (icav == kcav) {                        // in the same cavity modes
                                if (ilev == 0 || klev == 0) continue;  // skip S0
                                if (ilev == klev) {
                                    dV[idxdV] = KCI[j * Fmolm1 + ilev - 1];
                                } else if (j == 0) {
                                    dV[idxdV] = lcoeff;  // ci-modes
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int j = 0, idxdV = 0; j < N; ++j) {
            for (int i = 0; i < FF; ++i, ++idxdV) { dV[idxdV] *= sqrt(mod_W[j]); }
        }
    }
    first_call = false;

    if (flag < 2) return 0;

    return 0;
}

int PyrCav_ForceField::ForceField_epes_PC3(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    /**
     *
     * R[0]     : cavity-modes
     * R[1]     : molecule-modes
     * R[2...]  : molecule-modes
     *
     */

    int Fmolm1 = Fmol - 1;  // Fmolm1 == 2

    if (first_call)
        for (int i = 0; i < FF; ++i) V[i] = 0.0f;

    for (int i = 0, idx = 0; i < Fmol; ++i) {
        for (int k = 0; k < Fmol; ++k, ++idx) {
            if (i == 0 || k == 0) continue;  // skip S0
            if (i == k) {
                V[idx] = ECI[i];  // energies
                for (int j = 1, jk = i - 1; j < N; ++j, jk += Fmolm1) { V[idx] += KCI[jk] * sqrt(mod_W[j]) * R[j]; }
            } else {
                V[idx] = (                            // off-diagonal in S1-S2
                    gcoup * sqrt(mod_W[0]) * R[0]     // cavity-modes
                    + lcoeff * sqrt(mod_W[1]) * R[1]  // ci-modes
                );
            }
        }
    }

    if (flag < 1) return 0;

    if (first_call) {
        for (int i = 0; i < NFF; ++i) dV[i] = 0.0f;  // init to zero

        for (int j = 0, idxk = 0, idxdV = 0; j < N; ++j) {
            for (int i = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++idxdV) {
                    if (i == 0 || k == 0) continue;  // skip S0
                    if (j == 0) {
                        if (i != k) dV[idxdV] = gcoup;
                        continue;  // skip cavity-modes for KCI
                    }
                    if (i == k) {
                        dV[idxdV] = KCI[idxk++];
                    } else {
                        if (j == 1) dV[idxdV] = lcoeff;
                    }
                }
            }
        }
        for (int j = 0, idxdV = 0; j < N; ++j) {
            for (int i = 0; i < FF; ++i, ++idxdV) { dV[idxdV] *= sqrt(mod_W[j]); }
        }
    }
    first_call = false;

    if (flag < 2) return 0;

    return 0;
}

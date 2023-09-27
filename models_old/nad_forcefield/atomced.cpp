#include "atomced.h"

#include "../../utils/definitions.h"

AtomCED_ForceField::AtomCED_ForceField(const Param& iparm, const int& child)
    : Nad_ForceField(iparm){};  // for child to call
AtomCED_ForceField::AtomCED_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    Nb = Param_GetT(int, parm, "Nb");  // Num. for Discretization of one bath
    CHECK_EQ(Nb, N);                   // Nb*nbath must be N

    try {
        Param_GetV(ffflag, parm);
        fftype = AtomCED::_dict.at(ffflag);
    } catch (std::out_of_range& e) { LOG(FATAL) << e.what(); }


    // 1) initial system Hamiltonian
    ALLOCATE_PTR_TO_VECTOR(Hsys, FF);
    switch (fftype) {
        case AtomCED::CED2:
            CHECK_EQ(F, 2);
            Hsys[0] = -0.6738f;  // e1
            Hsys[1] = +1.034f;   // mu12
            Hsys[2] = +1.034f;   // mu21
            Hsys[3] = -0.2798f;  // e2
            Lcav    = 2.362e5;   // a.u.
            Rcav    = Lcav / 2;
            break;
        case AtomCED::CED3:
            CHECK_EQ(F, 3);
            Hsys[0] = -0.6738f;  // e1
            Hsys[1] = +1.034f;   // mu12
            Hsys[2] = +0.000f;   // mu13
            Hsys[3] = +1.034f;   // mu21
            Hsys[4] = -0.2798f;  // e2
            Hsys[5] = -2.536f;   // mu23
            Hsys[6] = +0.000f;   // mu31
            Hsys[7] = -2.536f;   // mu32
            Hsys[8] = -0.1547f;  // e3
            Lcav    = 2.362e5;   // a.u.
            Rcav    = Lcav / 2;
            break;
        default: {
            std::ifstream ifs(ffflag);
            num_real tmp;
            for (int i = 0; i < FF; ++i)
                if (ifs >> tmp) Hsys[i] = tmp / iou.ener;
            ifs.close();
            LOG(FATAL);
        }
    }

    // 2): initialize bath spectrum
    ALLOCATE_PTR_TO_VECTOR(omegas, Nb);
    ALLOCATE_PTR_TO_VECTOR(coeffs, Nb);
    for (int j = 0; j < Nb; ++j) {
        omegas[j] = (2 * j + 1) * lightspeed * phys::math::pi / Lcav;
        coeffs[j] = sqrt(2.0f / (epsilon0 * Lcav)) * sin((2 * j + 1) * phys::math::pi * Rcav / Lcav);
    }

    for (int j = 0; j < N; ++j) {
        mod_M[j] = 1.0f;
        mod_W[j] = omegas[j];

        double Qoverbeta;
        Qoverbeta     = 0.5f * mod_W[j];  // zero temp's Q
        mod_sigmaR[j] = std::sqrt(Qoverbeta / (mod_M[j] * mod_W[j] * mod_W[j]));
        mod_sigmaP[j] = std::sqrt(Qoverbeta * mod_M[j]);
    }

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

int AtomCED_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                        int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int AtomCED_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
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

int AtomCED_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                        const int& fdim) {
    // V
    for (int i = 0, idxV = 0; i < F; ++i)
        for (int k = 0; k < F; ++k, ++idxV) V[idxV] = (i == k) ? Hsys[idxV] : 0.0f;
    for (int j = 0; j < N; ++j) {
        for (int i = 0, idxV = 0; i < F; ++i) {
            for (int k = 0; k < F; ++k, ++idxV) {
                if (i == k) continue;
                V[idxV] += omegas[j] * R[j] * coeffs[j] * Hsys[idxV];
            }
        }
    }
    if (flag < 1) return 0;

    if (first_call) {
        // dV
        for (int j = 0, idxdV = 0; j < N; ++j) {
            for (int i = 0, idxV = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++idxV, ++idxdV) {
                    dV[idxdV] = (i == k) ? 0.0f : omegas[j] * coeffs[j] * Hsys[idxV];
                }
            }
        }
    }
    first_call = false;

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int AtomCED_ForceField::reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) {
    num_real* pdHj = dH;
    int Fadd1      = F + 1;
    double val     = 0.0f;
    for (int i = 0, idxd = 0; i < F; ++i, idxd += Fadd1) {
        for (int k = i + 1, idx1 = idxd + 1, idx2 = idxd + F; k < F; ++k, ++idx1, idx2 += F) {
            val += REAL_OF(rho[idx1]) * pdHj[idx2] + REAL_OF(rho[idx2]) * pdHj[idx1];
        }
    }
    double refval = coeffs[0] * omegas[0];
    for (int j = 0; j < N; ++j) fx[j] = val * coeffs[j] * omegas[j] / refval;
    return 0;
}

int AtomCED_ForceField::get_Nb() { return Nb; }

#include "liquidne_model.h"

#include "../../utils/elements.h"

extern "C" {
void liquidne_init_ccall(double*, const int&);
void liquidne_force_ccall(double*, double*, double*);
void liquidne_finalize_ccall();
}

//    void liquidne_init_ccall_none(double *cor, const int& p) {};
//    void liquidne_force_ccall_none(double *cor, double *V, double *dV) {};
//    void liquidne_finalize_ccall_none() {};

LiquidNe_ForceField::LiquidNe_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    CHECK_EQ(N, 648);

    const double m_Ne = ELEMENTS_MASS_NOAVG[10] / phys::au_2_amu;

    Pnbd = Param_GetT(int, parm, "nbead", 1);
    PtN  = Pnbd * N;

    for (int i = 0; i < N; ++i) mod_M[i] = m_Ne;
    for (int i = 0; i < N; ++i) { mod_R0[i] = 0.0f, mod_P0[i] = 0.0f, mod_sigmaR[i] = 0.0f, mod_sigmaP[i] = 0.0f; }

    tag = name() + "_" + tag;
}

LiquidNe_ForceField::~LiquidNe_ForceField() { liquidne_finalize_ccall(); };

int LiquidNe_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                         const int& rdim) {
    for (int i = 0; i < rdim; ++i) {
        // LOG(INFO) << i << ' ' << PtN;
        R[i] *= phys::au_2_ang * 1e-10L;
    }
    liquidne_force_ccall(R, V, dV);

    for (int i = 0; i < rdim; ++i) {
        R[i] /= phys::au_2_ang * 1e-10L;
        dV[i] = -dV[i] * phys::au_2_ang * 1e-10L;
        // dV[i] = 0.0;
        // std::cout << i << ' ' << rdim << ' ' << dV[i] << std::endl;
    }
    //
    return 0;
}

int LiquidNe_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    // BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    LOG(INFO) << "calling liquid init";
    liquidne_init_ccall(nr, rdim);
    // rand_gaussian(np, rdim);

    for (int i = 0; i < N; ++i) nm[i] = mod_M[i];
    for (int i = 0; i < rdim; ++i) nr[i] /= phys::au_2_ang * 1e-10L;
    for (int i = 0; i < rdim; ++i) np[i] = 0.0;
    return 0;
}

int LiquidNe_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return BO_ForceField ::ForceField_spec(nr, np, nm, rdim);
}

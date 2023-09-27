#include "spectrum_nadmodels.h"

#include "../forcefieldbase.h"
#include "lvcm_model.h"
#include "systembath.h"

Spectrum_NAD_ForceField::Spectrum_NAD_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    Param newparam = parm;
    Fminus1        = F - 1;

    Param_GetV(ground_shift, parm, 0.0f);
    ground_shift /= phys::au_2_ev;

    // build intrinsic Nad_ForceField
    newparam["fdim"] = (int) newparam["fdim"] - 1;
    pFF_inner        = new LVCM_ForceField(newparam);  // @note: global! don't locally declare!

    for (int j = 0; j < N; ++j) {  // copy nuclear DOFs
        mod_M[j]      = pFF_inner->mod_M[j];
        mod_W[j]      = pFF_inner->mod_W[j];
        mod_R0[j]     = pFF_inner->mod_R0[j];
        mod_P0[j]     = pFF_inner->mod_P0[j];
        mod_sigmaR[j] = pFF_inner->mod_sigmaR[j];
        mod_sigmaP[j] = pFF_inner->mod_sigmaP[j];
    }
    ALLOCATE_PTR_TO_VECTOR(workr_v, FF);
    ALLOCATE_PTR_TO_VECTOR(workr_dv, NFF);

    tag = name() + "_" + tag;
}

Spectrum_NAD_ForceField::~Spectrum_NAD_ForceField() { delete pFF_inner; }

int Spectrum_NAD_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                             int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int Spectrum_NAD_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                             const int& rdim) {
    pFF_inner->ForceField_npes(V, dV, ddV, R, P, flag, rdim);
    return 0;
}

int Spectrum_NAD_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag,
                                             const int& rdim, const int& fdim) {
    pFF_inner->ForceField_epes(V, dV, ddV, R, flag, rdim, Fminus1);
    // copy V from Fminus1-block to F-block
    for (int i = 0, idx1 = 0, idx2 = 0; i < F; ++i) {
        for (int k = 0; k < F; ++k, ++idx1) {
            if (i < Fminus1 && k < Fminus1) {
                workr_v[idx1] = V[idx2++];
            } else {
                workr_v[idx1] = 0.0f;
            }
        }
    }
    for (int i = 0; i < FF; ++i) V[i] = workr_v[i];
    V[FF - 1] = ground_shift;

    if (flag < 1) return 0;

    // copy dV from Fminus1-block to F-block
    if (first_call) {
        for (int j = 0, idx1 = 0, idx2 = 0; j < N; ++j) {
            for (int i = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++idx1) {
                    if (i < Fminus1 && k < Fminus1) {
                        workr_dv[idx1] = dV[idx2++];
                    } else {
                        workr_dv[idx1] = 0.0f;
                    }
                }
            }
        }
        for (int i = 0; i < NFF; ++i) dV[i] = workr_dv[i];
    }
    first_call = false;

    if (flag < 2) return 0;

    LOG(FATAL);

    return 0;
}

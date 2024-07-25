#include "ZnPc.h"

#include "../../utils/la_utils.h"
#include "../forcefieldbase.h"
#include "hamiltonian_data.h"
#include "systembath.h"


/**
 * @profile
 * about 55% time in EigenSolve problems & 20% time in update_rho
 */

using namespace ARRAY_EG;

ZnPc_ForceField::ZnPc_ForceField(const Param& iparm) : SystemBath_ForceField(iparm, 0) {
    Nb = sizeof(ZnPc_spectrum_data) / sizeof(ZnPc_spectrum_data[0]) / 5;
    CHECK_EQ(Nb, 165);
    CHECK_EQ(Nb, Param_GetT(int, parm, "Nb"));  // Num. for Discretization bath

    nmol  = Param_GetV(nmol, parm, 4);  // Num. for aggegated molecule
    nexc  = Param_GetV(nexc, parm, 2);  // Num. for excited state
    nbath = nmol;
    CHECK_EQ(nbath, Param_GetV(nbath, parm, 2));
    CHECK_EQ(F, nexc * nmol * nmol);  // nexc*nmol*nmol must be F
    CHECK_EQ(Nb * nbath, N);          // Nb*nmol must be N

    ALLOCATE_PTR_TO_VECTOR(idxarr_L, F);
    ALLOCATE_PTR_TO_VECTOR(idxarr_es, F);
    ALLOCATE_PTR_TO_VECTOR(idxarr_hs, F);
    ALLOCATE_PTR_TO_VECTOR(Hsys, FF);
    ALLOCATE_PTR_TO_VECTOR(eigen_T, FF);
    ALLOCATE_PTR_TO_VECTOR(eigen_E, F);
    ALLOCATE_PTR_TO_VECTOR(Etilde, nexc * nmol * nmol);
    ALLOCATE_PTR_TO_VECTOR(Vtilde, nexc * nmol * nmol);
    ALLOCATE_PTR_TO_VECTOR(te_tilde, nexc * nmol * nmol);
    ALLOCATE_PTR_TO_VECTOR(th_tilde, nexc * nmol * nmol);
    ALLOCATE_PTR_TO_VECTOR(tect_tilde, nexc * nmol * nmol);
    ALLOCATE_PTR_TO_VECTOR(thct_tilde, nexc * nmol * nmol);

    init_Hamiltonian();  // build Hsys and idxarrs

    ALLOCATE_PTR_TO_VECTOR(omegas, N);
    ALLOCATE_PTR_TO_VECTOR(coeffs, N);
    ALLOCATE_PTR_TO_VECTOR(dQe1, Nb);
    ALLOCATE_PTR_TO_VECTOR(dQe2, Nb);
    ALLOCATE_PTR_TO_VECTOR(dQc, Nb);
    ALLOCATE_PTR_TO_VECTOR(dQa, Nb);
    ALLOCATE_PTR_TO_VECTOR(w2dQe1, Nb);
    ALLOCATE_PTR_TO_VECTOR(w2dQe2, Nb);
    ALLOCATE_PTR_TO_VECTOR(w2dQc, Nb);
    ALLOCATE_PTR_TO_VECTOR(w2dQa, Nb);

    ALLOCATE_PTR_TO_VECTOR(Xnj, NFF);

    // for fast sse
    Param_Reset(L, 4);  // 4-nonzero values in each Q
    ALLOCATE_PTR_TO_VECTOR(QL, L * nbath * FF);
    ALLOCATE_PTR_TO_VECTOR(CL, L * Nb);

    for (int i = 0, idx = 0; i < Nb; ++i) {
        omegas[i] = ZnPc_spectrum_data[idx++] / phys::au_2_wn;
        dQe1[i]   = ZnPc_spectrum_data[idx++];
        dQe2[i]   = ZnPc_spectrum_data[idx++];
        dQc[i]    = ZnPc_spectrum_data[idx++];
        dQa[i]    = ZnPc_spectrum_data[idx++];
        // @OPT
        w2dQe1[i]      = omegas[i] * omegas[i] * dQe1[i];
        w2dQe2[i]      = omegas[i] * omegas[i] * dQe2[i];
        w2dQa[i]       = omegas[i] * omegas[i] * dQa[i];
        w2dQc[i]       = omegas[i] * omegas[i] * dQc[i];
        CL[0 * Nb + i] = w2dQe1[i];
        CL[1 * Nb + i] = w2dQe2[i];
        CL[2 * Nb + i] = w2dQa[i];
        CL[3 * Nb + i] = w2dQc[i];
    }
    // LOG(FATAL);
    for (int ibath = 0, idx = 0; ibath < nbath; ++ibath) {
        for (int j = 0; j < Nb; ++j, ++idx) {
            mod_M[idx] = 1.0f, mod_W[idx] = omegas[j];
            for (int i = 0; i < F; ++i) {
                int iL = idxarr_L[i], ie = idxarr_es[i], ih = idxarr_hs[i];
                double& Xval = Xnj[idx * FF + i * (F + 1)];
                if (ibath == ie && ibath == ih) {
                    Xval   = (iL == 0) ? w2dQe1[j] : w2dQe2[j];
                    int ii = (iL == 0) ? 0 : 1;
                    // QL
                    QL[ii * nbath * FF + ibath * FF + i * (F + 1)] = 1;
                } else if (ibath == ie && ibath != ih) {
                    Xval = w2dQa[j];
                    // QL
                    QL[2 * nbath * FF + ibath * FF + i * (F + 1)] = 1;
                } else if (ibath != ie && ibath == ih) {
                    Xval = w2dQc[j];
                    // QL
                    QL[3 * nbath * FF + ibath * FF + i * (F + 1)] = 1;
                } else {
                    Xval = 0.0f;
                }
            }
        }
    }

    // init bath
    mybath = new Bath(parm, nbath, F);

    BO_ForceField::ForceField_init_default_build(mybath->beta, N);
    CheckForceField();
    tag = name() + "_" + tag;
}

int ZnPc_ForceField::index(const int& iL, const int& ie, const int& ih) {
    int pos;
    if (ie == ih) { pos = iL * nmol * nmol + ie; }
    if (ie < ih) {
        int l = ih - ie;
        pos   = iL * nmol * nmol + nmol + (2 * nmol - l) * (l - 1) / 2 + ie;
    }
    if (ih < ie) {
        int l = ie - ih;
        pos   = iL * nmol * nmol + nmol + (nmol - 1) * nmol / 2 + (2 * nmol - l) * (l - 1) / 2 + ih;
    }
    return pos;
}

int ZnPc_ForceField::init_Hamiltonian() {
    // initial to zero
    for (int i = 0; i < F; ++i)
        idxarr_L[i] = 0, idxarr_es[i] = 0, idxarr_hs[i] = 0, eigen_E[i] = 0.0f;  // excited number
    for (int i = 0; i < FF; ++i) Hsys[i] = 0.0f, eigen_T[i] = 0.0f;
    for (int i = 0; i < nexc * nmol * nmol; ++i) Etilde[i] = 0.0f;
    for (int i = 0; i < nexc * nexc * nmol; ++i) {
        Vtilde[i] = 0, te_tilde[i] = 0, th_tilde[i] = 0, tect_tilde[i] = 0, thct_tilde[i] = 0;
    }

    // record index & calc Etilde
    for (int iL = 0; iL < nexc; ++iL) {
        for (int ie = 0; ie < nmol; ++ie) {
            for (int ih = 0; ih < nmol; ++ih) {
                int i       = index(iL, ie, ih);
                idxarr_L[i] = iL, idxarr_es[i] = ie, idxarr_hs[i] = ih;

                int l      = std::abs(ie - ih);
                double& Ex = Etilde[iL * nmol * nmol + ie * nmol + ih];
                if (l == 0) {  // excited state
                    Ex = (iL == 0) ? 1.88f : 1.94f;
                    if (ie == 0 || ie == nmol - 1) Ex = (iL == 0) ? 1.90f : 1.96f;
                } else if (l == 1) {  // CT state
                    Ex = (iL == 0) ? 2.30f : 2.39f;
                    if (ih == 0 || ih == nmol - 1) Ex = (iL == 0) ? 2.44f : 2.53f;
                    if (ie == 0 || ie == nmol - 1) Ex = (iL == 0) ? 2.18f : 2.27f;
                } else if (l == 2) {                 // CS state
                    Ex = (iL == 0) ? 3.20f : 3.28f;  //@ERROR
                    if (ih == 0 || ih == nmol - 1) Ex = (iL == 0) ? 3.33f : 3.41f;
                    if (ie == 0 || ie == nmol - 1) Ex = (iL == 0) ? 3.08f : 3.15f;
                } else if (l == 3) {
                    Ex = (iL == 0) ? 3.67f : 3.74f;
                } else {
                    Ex = (iL == 0) ? -1.930f / l + 4.228f : -1.900f / l + 4.288f;
                }
            }
        }
    }
    for (int i = 0; i < F; ++i) {
        LOG(WARNING) << "idx: " << i << " (L, e, h) = (" << idxarr_L[i] << ", " << idxarr_es[i] << ", " << idxarr_hs[i]
                     << ")" << std::endl;
    }
    int D0 = nexc * nmol, D1 = nmol;

    Vtilde[1] = 0.184f, Vtilde[2] = 0.0508f, Vtilde[3] = 0.020f;  // V11
    Vtilde[D0 + D1 + 1] = 0.169f, Vtilde[D0 + D1 + 2] = 0.0476f,
                     Vtilde[D0 + D1 + 3] = 0.0189f;  // V22
    // V12 = V21 = 0

    te_tilde[1] = -0.0172f, te_tilde[2] = -0.0176f;  // @ERROR -0.0177
    te_tilde[D0 + D1 + 1] = 0.0269f, te_tilde[D0 + D1 + 2] = -0.0122f;
    te_tilde[D1 + 1] = -0.136f;  //@ERROR -0.1363
    te_tilde[D0 + 1] = -0.133f;

    th_tilde[1] = -0.0126f, th_tilde[2] = 0.022f;
    th_tilde[D0 + D1 + 1] = -0.0160f, th_tilde[D0 + D1 + 2] = 0.0216f;
    th_tilde[D0 + 1] = 0.0140f;  // only off-diagonal (th_[L2L1])

    tect_tilde[1] = -0.0144f, tect_tilde[2] = -0.014f;
    tect_tilde[D0 + D1 + 1] = -0.0f, tect_tilde[D0 + D1 + 2] = -0.0108f;
    tect_tilde[D1 + 1] = -0.133f;  // @ERROR -0.1327
    tect_tilde[D0 + 1] = -0.133f;  // @ERROR -0.1327

    thct_tilde[1] = -0.025f, thct_tilde[2] = 0.0196;
    thct_tilde[D0 + D1 + 1] = -0.025f, thct_tilde[D0 + D1 + 2] = 0.0195f;
    // none off-diagonal

    for (int i1 = 0; i1 < F; ++i1) {
        int iL1 = idxarr_L[i1], ie1 = idxarr_es[i1], ih1 = idxarr_hs[i1];
        int l1 = std::abs(ie1 - ih1);
        for (int i2 = 0; i2 < F; ++i2) {
            int iL2 = idxarr_L[i2], ie2 = idxarr_es[i2], ih2 = idxarr_hs[i2];
            int l2 = std::abs(ie2 - ih2);

            num_real& Hx = Hsys[i1 * F + i2];
            // state energy
            if (i1 == i2) {
                Hx = Etilde[iL1 * nmol * nmol + ie1 * nmol + ih1];
                // couplings
            } else if (ie1 == ih1 && ie2 == ih2 && ie1 != ie2) {
                Hx = Vtilde[iL1 * nexc * nmol + iL2 * nmol + std::abs(ie1 - ie2)];
            } else if (ie1 == ih1 && ie2 != ih2 && ih1 == ih2) {
                Hx = te_tilde[iL1 * nexc * nmol + iL2 * nmol + std::abs(ie1 - ie2)];
            } else if (ie1 != ih1 && ie2 == ih2 && ih1 == ih2) {
                Hx = te_tilde[iL2 * nexc * nmol + iL1 * nmol + std::abs(ie1 - ie2)];
            } else if (ie1 == ih1 && ie2 != ih2 && ie1 == ie2) {
                Hx = th_tilde[iL1 * nexc * nmol + iL2 * nmol + std::abs(ih1 - ih2)];
            } else if (ie1 != ih1 && ie2 == ih2 && ie1 == ie2) {
                Hx = th_tilde[iL2 * nexc * nmol + iL1 * nmol + std::abs(ih1 - ih2)];
            } else if (ie1 != ih1 && ie2 != ih2 && ih1 == ih2 && ie1 < ie2) {
                Hx = tect_tilde[iL1 * nexc * nmol + iL2 * nmol + std::abs(ie1 - ie2)];
            } else if (ie1 != ih1 && ie2 != ih2 && ih1 == ih2 && ie1 > ie2) {
                Hx = tect_tilde[iL2 * nexc * nmol + iL1 * nmol + std::abs(ie1 - ie2)];
            } else if (ie1 != ih1 && ie2 != ih2 && ie1 == ie2 && ih1 < ih2) {
                Hx = thct_tilde[iL1 * nexc * nmol + iL2 * nmol + std::abs(ih1 - ih2)];
            } else if (ie1 != ih1 && ie2 != ih2 && ie1 == ie2 && ih1 > ih2) {
                Hx = thct_tilde[iL2 * nexc * nmol + iL1 * nmol + std::abs(ih1 - ih2)];
            }
        }
    }
    for (int i = 0; i < FF; ++i) Hsys[i] /= phys::au_2_ev;
    return 0;
}

int ZnPc_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                     int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int ZnPc_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                     const int& rdim) {
    plFunction();

    double Vnucl = 0.0f;
    for (int j = 0; j < N; ++j) {  // @profiling: self-time from 20% => 13%
        double mwRj = mod_M[j] * mod_W[j] * R[j];
        Vnucl += (P[j] * P[j] + mwRj * mwRj) / mod_M[j];
        dV[j] = mod_W[j] * mwRj;
    }
    V[0] = 0.5f * Vnucl;

    if (flag < 2) return 0;
    ARRAY_CLEAR(ddV, NN);
    for (int j = 0, idx = 0, add = N + 1; j < N; ++j, idx += add) ddV[idx] = mod_M[j] * mod_W[j] * mod_W[j];
    return 0;
}

int ZnPc_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                     const int& fdim, const int& itraj, const int& isamp) {
    plFunction();
    int iL, ie, ih;

    // V: about 50% faster than above [Nb*F] order less than [N*F] order
    for (int i = 0; i < FF; ++i) V[i] = Hsys[i];
    int Fadd1 = F + 1;
    for (int i = 0, idxV = 0; i < F; ++i, idxV += Fadd1) {
        iL = idxarr_L[i], ie = idxarr_es[i], ih = idxarr_hs[i];
        if (ie == ih) {
            for (int j = 0, idxR = ie * Nb; j < Nb; ++j, ++idxR) {
                V[idxV] -= ((iL == 0) ? w2dQe1[j] : w2dQe2[j]) * R[idxR];
            }
        } else {
            for (int j = 0, idxRe = ie * Nb, idxRh = ih * Nb; j < Nb; ++j, ++idxRe, ++idxRh) {
                V[idxV] -= (w2dQa[j] * R[idxRe] + w2dQc[j] * R[idxRh]);
            }
        }
    }
    if (flag < 1) return 0;

    // dV: similar time ? [Nb*F] order, only calc once time
    if (isamp < 1) {
        for (int i = 0; i < F; ++i) {
            iL = idxarr_L[i], ie = idxarr_es[i], ih = idxarr_hs[i];
            if (ie == ih) {
                int idxdV = ie * Nb * FF + i * (F + 1);
                for (int j = 0; j < Nb; ++j, idxdV += FF) { dV[idxdV] = -((iL == 0) ? w2dQe1[j] : w2dQe2[j]); }
            } else {
                int idxdVe = ie * Nb * FF + i * (F + 1), idxdVh = ih * Nb * FF + i * (F + 1);
                for (int j = 0; j < Nb; ++j, idxdVe += FF, idxdVh += FF) {
                    dV[idxdVe] = -w2dQa[j];
                    dV[idxdVh] = -w2dQc[j];
                }
            }
        }
    }
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int ZnPc_ForceField::get_nbath() { return nbath; }
int ZnPc_ForceField::get_Nb() { return Nb; }

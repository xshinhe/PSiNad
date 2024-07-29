#include "kids/Kernel_Update_p.h"

#include "kids/Kernel_Monodromy.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_p::getName() { return "Kernel_Update_p"; }

int Kernel_Update_p::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_p::setInputParam_impl(std::shared_ptr<Param> PM) {
    use_smooth = _param->get_bool({"solver.use_smooth"}, LOC(), false);
}

void Kernel_Update_p::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    f       = DS->def(DATA::integrator::f);
    fadd    = DS->def(DATA::integrator::fadd);
    p       = DS->def(DATA::integrator::p);
    ve      = DS->def(DATA::integrator::ve);
    minv    = DS->def(DATA::integrator::minv);
    mono    = DS->def(DATA::integrator::monodromy::mono);
    monodt  = DS->def(DATA::integrator::monodromy::monodt);
    MFFtmp1 = DS->def(DATA::integrator::monodromy::MFFtmp1);
    MFFtmp2 = DS->def(DATA::integrator::monodromy::MFFtmp2);
    MFFtmp3 = DS->def(DATA::integrator::monodromy::MFFtmp3);
    MFFtmp4 = DS->def(DATA::integrator::monodromy::MFFtmp4);
    MFFtmp5 = DS->def(DATA::integrator::monodromy::MFFtmp5);
    MFFtmp6 = DS->def(DATA::integrator::monodromy::MFFtmp6);
    mask    = DS->def(DATA::integrator::forceeval::mask);
    dmask   = DS->def(DATA::integrator::forceeval::dmask);
    T       = DS->def(DATA::model::rep::T);
    grad    = DS->def(DATA::model::grad);
    hess    = DS->def(DATA::model::hess);
    dV      = DS->def(DATA::model::dV);
    ddV     = DS->def(DATA::model::ddV);
    dE      = DS->def(DATA::model::rep::dE);
    nac     = DS->def(DATA::model::rep::nac);
    eig     = DS->def(DATA::model::rep::eig);
    rho_nuc = DS->def(DATA::integrator::rho_nuc);
    c       = DS->def(DATA::integrator::c);
    Ekin    = DS->def(DATA::integrator::Ekin);
    dt      = DS->def(DATA::flowcontrol::dt);
}

Status& Kernel_Update_p::initializeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;
        Ekin[0]         = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
    }
    return stat;
}

Status& Kernel_Update_p::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    if (!use_smooth) {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            kids_real* f = this->f + iP * Dimension::N;
            kids_real* p = this->p + iP * Dimension::N;
            for (int i = 0; i < Dimension::N; ++i) { p[i] -= f[i] * scale * dt[0]; }
        }
    } else {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            kids_real*    f       = this->f + iP * Dimension::N;
            kids_real*    p       = this->p + iP * Dimension::N;
            kids_real*    T       = this->T + iP * Dimension::FF;        //
            kids_real*    dE      = this->dE + iP * Dimension::NFF;      //
            kids_real*    eig     = this->eig + iP * Dimension::F;       //
            kids_complex* mask    = this->mask + iP * Dimension::FF;     //
            kids_complex* dmask   = this->dmask + iP * Dimension::NFF;   //
            kids_complex* rho_nuc = this->rho_nuc + iP * Dimension::FF;  //

            kids_complex im = phys::math::im;
            for (int i = 0, ii = 0, ik = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1, ++ik) {
                    kids_real    dlam    = (eig[i] - eig[k]);
                    kids_complex idlamdt = im * dlam * dt[0];
                    kids_complex expterm = std::exp(idlamdt);
                    mask[ik]             = (i == k) ? 1.0e0 : (expterm - 1.0e0) / idlamdt;
                    for (int J = 0, Jik = ik, Jii = ii, Jkk = kk; J < Dimension::N;
                         ++J, Jik += Dimension::FF, Jii += Dimension::FF, Jkk += Dimension::FF) {
                        dmask[Jik] = (i == k)
                                         ? 0.0e0
                                         : (dlam * dt[0] * expterm + im * (expterm - 1.0e0)) / (dlam * dlam * dt[0]);
                        dmask[Jik] *= (dE[Jii] - dE[Jkk]);
                    }
                }
            }
            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::nuc_repr_type,  //
                                             SpacePolicy::L);
            for (int I = 0; I < Dimension::N; ++I) {
                kids_real* dEI = dE + I * Dimension::FF;
                p[I] -= fadd[I] * scale * dt[0];  //
                p[I] -= grad[I] * scale * dt[0];  //
                for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp1[ik] = mask[ik] * dEI[ik];
                ARRAY_MATMUL3_TRANS2(MFFtmp1, T, MFFtmp1, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                p[I] -= std::real(ARRAY_TRACE2(rho_nuc, MFFtmp1, Dimension::F, Dimension::F)) * scale * dt[0];
            }
            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::nuc_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
        }
    }

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;
        Ekin[0]         = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) { Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i]; }
    }
    // trace on monodromy
    if (Kernel_Monodromy::enable) update_monodromy();
    return stat;
}

void Kernel_Update_p::update_monodromy() {
    int N0   = 0;
    int N1   = Dimension::N;
    int N2   = Dimension::N + Dimension::F;
    int N3   = 2 * Dimension::N + Dimension::F;
    int N4   = 2 * Dimension::N + 2 * Dimension::F;
    int N4N4 = N4 * N4;
    if (!use_smooth) {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            kids_real*    mono    = this->mono + iP * N4N4;              //
            kids_real*    monodt  = this->monodt + iP * N4N4;            //
            kids_real*    hess    = this->hess + iP * Dimension::NN;     //
            kids_real*    dV      = this->dV + iP * Dimension::NFF;      //
            kids_real*    ddV     = this->ddV + iP * Dimension::NNFF;    //
            kids_complex* rho_nuc = this->rho_nuc + iP * Dimension::FF;  //
            kids_complex* c       = this->c + iP * Dimension::F;         // @bug for rep?

            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::nuc_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_EYE(monodt, N4);
            for (int I = 0, IJ = 0; I < Dimension::N; ++I) {
                kids_real* dVI = dV + I * Dimension::FF;
                ARRAY_MATMUL(MFFtmp2, dVI, c, Dimension::F, Dimension::F, 1);
                for (int k = 0; k < Dimension::F; ++k) {
                    monodt[(N2 + I) * N4 + (N1 + k)] = 2.0e0 * std::real(MFFtmp2[k]) * scale * dt[0];
                    monodt[(N2 + I) * N4 + (N3 + k)] = 2.0e0 * std::imag(MFFtmp2[k]) * scale * dt[0];
                }
                for (int J = 0; J < Dimension::N; ++J, ++IJ) {  //
                    kids_real* ddVIJ                 = ddV + IJ * Dimension::FF;
                    monodt[(N2 + I) * N4 + (N0 + J)] = -hess[IJ] * scale * dt[0];
                    monodt[(N2 + I) * N4 + (N0 + J)] +=
                        -std::real(ARRAY_TRACE2(ddVIJ, rho_nuc, Dimension::F, Dimension::F)) * scale * dt[0];
                }
            }
            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::nuc_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_MATMUL(mono, monodt, mono, N4, N4, N4);
        }
    } else {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            kids_real*    mono    = this->mono + iP * N4N4;              //
            kids_real*    monodt  = this->monodt + iP * N4N4;            //
            kids_real*    hess    = this->hess + iP * Dimension::NN;     //
            kids_real*    dV      = this->dV + iP * Dimension::NFF;      //
            kids_real*    T       = this->T + iP * Dimension::FF;        //
            kids_real*    eig     = this->eig + iP * Dimension::F;       //
            kids_real*    dE      = this->dE + iP * Dimension::NFF;      //
            kids_real*    nac     = this->nac + iP * Dimension::NFF;     //
            kids_real*    ddV     = this->ddV + iP * Dimension::NNFF;    //
            kids_complex* mask    = this->mask + iP * Dimension::FF;     //
            kids_complex* dmask   = this->dmask + iP * Dimension::NFF;   //
            kids_complex* c       = this->c + iP * Dimension::F;         // @bug for rep?
            kids_complex* rho_nuc = this->rho_nuc + iP * Dimension::FF;  //

            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::nuc_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_EYE(monodt, N4);
            for (int I = 0, IJ = 0; I < Dimension::N; ++I) {
                kids_real* dEI = dE + I * Dimension::FF;
                for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp1[ik] = mask[ik] * dEI[ik];
                ARRAY_MATMUL3_TRANS2(MFFtmp1, T, MFFtmp1, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                ARRAY_MATMUL(MFFtmp2, MFFtmp1, c, Dimension::F, Dimension::F, 1);
                for (int k = 0; k < Dimension::F; ++k) {
                    monodt[(N2 + I) * N4 + (N1 + k)] = 2.0e0 * std::real(MFFtmp2[k]) * scale * dt[0];
                    monodt[(N2 + I) * N4 + (N3 + k)] = 2.0e0 * std::imag(MFFtmp2[k]) * scale * dt[0];
                }

                for (int J = 0; J < Dimension::N; ++J, ++IJ) {  //
                    kids_real*    ddVIJ  = ddV + IJ * Dimension::FF;
                    kids_real*    nacJ   = nac + J * Dimension::FF;
                    kids_complex* dmaskJ = dmask + J * Dimension::FF;

                    monodt[(N2 + I) * N4 + (N0 + J)] = -hess[IJ] * scale * dt[0];

                    for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp1[ik] = mask[ik] * dEI[ik];
                    ARRAY_MATMUL(MFFtmp1, nacJ, MFFtmp1, Dimension::F, Dimension::F, Dimension::F);

                    for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp2[ik] = mask[ik] * dEI[ik];
                    ARRAY_MATMUL(MFFtmp2, MFFtmp2, nacJ, Dimension::F, Dimension::F, Dimension::F);

                    ARRAY_MATMUL(MFFtmp3, dEI, nacJ, Dimension::F, Dimension::F, Dimension::F);
                    for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp3[ik] *= mask[ik];

                    ARRAY_MATMUL(MFFtmp4, nacJ, dEI, Dimension::F, Dimension::F, Dimension::F);
                    for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp4[ik] *= mask[ik];

                    ARRAY_MATMUL3_TRANS1(MFFtmp5, T, ddVIJ, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                    for (int ik = 0; ik < Dimension::FF; ++ik) MFFtmp5[ik] *= mask[ik];

                    for (int ik = 0; ik < Dimension::FF; ++ik)
                        MFFtmp6[ik] =
                            dE[ik] * dmaskJ[ik] + MFFtmp5[ik] + MFFtmp1[ik] - MFFtmp2[ik] + MFFtmp3[ik] - MFFtmp4[ik];
                    ARRAY_MATMUL3_TRANS2(MFFtmp6, T, MFFtmp6, T, Dimension::F, Dimension::F, Dimension::F,
                                         Dimension::F);

                    monodt[(N2 + I) * N4 + (N0 + J)] +=
                        -std::real(ARRAY_TRACE2(MFFtmp6, rho_nuc, Dimension::F, Dimension::F)) * scale * dt[0];
                }
            }
            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::nuc_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_MATMUL(mono, monodt, mono, N4, N4, N4);
        }
    }
}

};  // namespace PROJECT_NS
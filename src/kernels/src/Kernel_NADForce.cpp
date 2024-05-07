#include "kids/Kernel_NADForce.h"

#include "kids/Kernel_Elec.h"
#include "kids/Kernel_Representation.h"
#include "kids/Kernel_Update.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_NADForce::getName() { return "Kernel_NADForce"; }

int Kernel_NADForce::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_NADForce::setInputParam_impl(std::shared_ptr<Param>& PM) {
    FORCE_OPT::BATH_FORCE_BILINEAR = _param->get_bool("BATH_FORCE_BILINEAR", LOC(), false);
    offd_projected                 = _param->get_bool("offd_projected", LOC(), true);
};

void Kernel_NADForce::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    f        = DS->def(DATA::integrator::f);
    p        = DS->def(DATA::integrator::p);
    m        = DS->def(DATA::integrator::m);
    fadd     = DS->def(DATA::integrator::fadd);
    fproj    = DS->def(DATA::integrator::tmp::fproj);
    grad     = DS->def(DATA::model::grad);
    dV       = DS->def(DATA::model::dV);
    dE       = DS->def(DATA::model::rep::dE);
    T        = DS->def(DATA::model::rep::T);
    succ_ptr = DS->def(DATA::iter::succ);

    switch (Kernel_Representation::nuc_repr_type) {
        case RepresentationPolicy::Diabatic:
            Force = dV;
            break;
        case RepresentationPolicy::Adiabatic:
            Force = dE;
            break;
    }
}

Status& Kernel_NADForce::initializeKernel_impl(Status& stat) { return executeKernel(stat); }

Status& Kernel_NADForce::executeKernel_impl(Status& stat) {
    if (!succ_ptr[0]) return stat;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        int*          occ_nuc = Kernel_Elec::occ_nuc + iP;
        kids_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;

        kids_real* f     = this->f + iP * Dimension::N;
        kids_real* p     = this->p + iP * Dimension::N;
        kids_real* m     = this->m + iP * Dimension::N;
        kids_real* fadd  = this->fadd + iP * Dimension::N;
        kids_real* grad  = this->grad + iP * Dimension::N;
        kids_real* Force = this->Force + iP * Dimension::NFF;
        kids_real* T     = this->T + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////
        switch (NADForce_type) {
            case NADForcePolicy::BO: {
                for (int j = 0, idxdV0 = 0; j < Dimension::N; ++j, idxdV0 += Dimension::FF)
                    f[j] = Force[j * Dimension::FF + (*occ_nuc) * Dimension::Fadd1];
                break;
            }
            case NADForcePolicy::EHR: {
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::inp_repr_type,
                                                 Kernel_Representation::nuc_repr_type, SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {  // for both dV & dE (only for FMO-like model)
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        double* Forceb0 = Force + b0FF;
                        double  fb0     = std::real(ARRAY_TRACE2(rho_nuc, Forceb0, Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj] = fb0 * Force[bjbb] / Force[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        double* dVj = Force + jFF;
                        f[j]        = std::real(ARRAY_TRACE2(rho_nuc, dVj, Dimension::F, Dimension::F));
                    }
                }
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::nuc_repr_type,
                                                 Kernel_Representation::inp_repr_type,
                                                 SpacePolicy::L);  // not need
                break;
            }
            case NADForcePolicy::CV: {
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::inp_repr_type,
                                                 Kernel_Representation::nuc_repr_type, SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {  // for both dV & dE (only for FMO-like model)
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        double* Forceb0 = Force + b0FF;
                        double  fb0     = Forceb0[(*occ_nuc) * Dimension::Fadd1];
                        double  fprojb0 = std::real(ARRAY_TRACE2_OFFD(rho_nuc, Forceb0, Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj]     = fb0 * Force[bjbb] / Force[b0bb];
                            fproj[bj] = fprojb0 * Force[bjbb] / Force[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        double* dVj = Force + jFF;
                        f[j]        = dVj[(*occ_nuc) * Dimension::Fadd1];
                        fproj[j]    = std::real(ARRAY_TRACE2_OFFD(rho_nuc, dVj, Dimension::F, Dimension::F));
                    }
                }
                if (offd_projected) {  // then the offdiagonal force is projected
                    double fdotR = 0.0e0, PdotR = 0.0e0;
                    for (int j = 0; j < Dimension::N; ++j) fdotR += fproj[j] * p[j] / m[j], PdotR += p[j] * p[j] / m[j];
                    for (int j = 0; j < Dimension::N; ++j) fproj[j] -= fdotR / PdotR * p[j];
                }
                for (int j = 0; j < Dimension::N; ++j) f[j] += fproj[j];

                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::nuc_repr_type,
                                                 Kernel_Representation::inp_repr_type,
                                                 SpacePolicy::L);  // not need
                break;
            }
            case NADForcePolicy::BOSD:
            case NADForcePolicy::CVSD: {  // smooth dynamics
                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::inp_repr_type,
                                                 Kernel_Representation::nuc_repr_type, SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {  // for both dV & dE (only for FMO-like model)
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        double* Forceb0 = Force + b0FF;
                        double  fb0     = std::real(ARRAY_TRACE2_DIAG(rho_nuc, Forceb0, Dimension::F, Dimension::F));
                        double  fprojb0 = std::real(ARRAY_TRACE2_OFFD(rho_nuc, Forceb0, Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj]     = fb0 * Force[bjbb] / Force[b0bb];
                            fproj[bj] = fprojb0 * Force[bjbb] / Force[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        double* dVj = Force + jFF;
                        f[j]        = std::real(ARRAY_TRACE2_DIAG(rho_nuc, dVj, Dimension::F, Dimension::F));
                        fproj[j]    = std::real(ARRAY_TRACE2_OFFD(rho_nuc, dVj, Dimension::F, Dimension::F));
                    }
                }

                if (NADForce_type == NADForcePolicy::CVSD) {
                    if (offd_projected) {
                        double fdotR = 0.0e0, vnorm2 = 0.0e0;
                        for (int j = 0; j < Dimension::N; ++j) fdotR += fproj[j] * p[j], vnorm2 += p[j] * p[j];
                        if (vnorm2 > 1.0e-8) {
                            for (int j = 0; j < Dimension::N; ++j) fproj[j] += (fproj[j] - fdotR / vnorm2 * p[j]);
                        }
                    } else {
                        for (int j = 0; j < Dimension::N; ++j) fproj[j] += fproj[j];
                    }
                }

                for (int j = 0; j < Dimension::N; ++j) f[j] += fproj[j];

                Kernel_Representation::transform(rho_nuc, T, Dimension::F, Kernel_Representation::nuc_repr_type,
                                                 Kernel_Representation::inp_repr_type,
                                                 SpacePolicy::L);  // not need
                break;
            }
        }
        for (int j = 0; j < Dimension::N; ++j) f[j] += grad[j];
        for (int j = 0; j < Dimension::N; ++j) f[j] += fadd[j];  // additional force
    }
    return stat;
}

NADForcePolicy::_type Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;

namespace FORCE_OPT {
int  nbath               = 1;
int  Nb                  = 1;
bool BATH_FORCE_BILINEAR = false;
};  // namespace FORCE_OPT

};  // namespace PROJECT_NS

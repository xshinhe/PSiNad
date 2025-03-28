#include "kids/Kernel_NAForce.h"

#include <algorithm>

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

double calc_alpha(kids_real* V, int i = 0, int k = 1, int F = 2) {  // acoording to mix angle
    int ii = i * (F + 1), kk = k * (F + 1), ik = i * F + k;
    if (V[ii] == V[kk]) return 1.0e0;
    double res = 2.0e0 / phys::math::pi * atan(2 * V[ik] / (V[ii] - V[kk]));
    if (abs(res) < 1.0e-4) res = copysign(1.0e-4, res);
    return res;
}

const std::string Kernel_NAForce::getName() { return "Kernel_NAForce"; }

int Kernel_NAForce::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_NAForce::setInputParam_impl(std::shared_ptr<Param> PM) {
    NAForce_type                   = NAForcePolicy::_dict.at(  //
        _param->get_string({"solver.naforce"}, LOC(), "EHR"));
    FORCE_OPT::BATH_FORCE_BILINEAR = _param->get_bool({"solver.BATH_FORCE_BILINEAR"}, LOC(), false);
    offd_projected                 = _param->get_bool({"solver.offd_projected"}, LOC(), true);
};

void Kernel_NAForce::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    dt_ptr = DS->def(DATA::flowcontrol::dt);
    f      = DS->def(DATA::integrator::f);
    // fmat     = DS->def(DATA::integrator::fmat);
    p        = DS->def(DATA::integrator::p);
    m        = DS->def(DATA::integrator::m);
    grad     = DS->def(DATA::model::grad);
    V        = DS->def(DATA::model::V);
    dV       = DS->def(DATA::model::dV);
    dE       = DS->def(DATA::model::rep::dE);
    T        = DS->def(DATA::model::rep::T);
    succ_ptr = DS->def(DATA::flowcontrol::succ);

    Epot = DS->def(DATA::integrator::Epot);
    vpes = DS->def(DATA::model::vpes);

    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    rho_ele = DS->def(DATA::integrator::rho_ele);
    rho_nuc = DS->def(DATA::integrator::rho_nuc);

    // Adaptive or Projective
    alpha = DS->def(DATA::integrator::alpha);
    fadd  = DS->def(DATA::integrator::fadd);
    fproj = DS->def(DATA::integrator::tmp::fproj);
    ftmp  = DS->def(DATA::integrator::tmp::ftmp);
    wrho  = DS->def(DATA::integrator::tmp::wrho);

    switch (Kernel_Representation::nuc_repr_type) {
        case RepresentationPolicy::Diabatic:
            EMat     = DS->def(DATA::model::V);
            ForceMat = DS->def(DATA::model::dV);
            break;
        case RepresentationPolicy::Adiabatic:
            EMat     = DS->def(DATA::model::rep::E);
            ForceMat = DS->def(DATA::model::rep::dE);
            break;
    }
}

Status& Kernel_NAForce::initializeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto V     = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        auto alpha = this->alpha.subspan(iP, 1);
        alpha[0]   = calc_alpha(V.data());
    }
    executeKernel(stat);
    return stat;
}

Status& Kernel_NAForce::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    for (int iP = 0; iP < Dimension::P_NOW; ++iP) {
        auto occ_nuc  = this->occ_nuc.subspan(iP, 1);
        auto rho_ele  = this->rho_ele.subspan(iP * Dimension::FF, Dimension::FF);
        auto rho_nuc  = this->rho_nuc.subspan(iP * Dimension::FF, Dimension::FF);
        auto f        = this->f.subspan(iP * Dimension::N, Dimension::N);
        auto p        = this->p.subspan(iP * Dimension::N, Dimension::N);
        auto m        = this->m.subspan(iP * Dimension::N, Dimension::N);
        auto fadd     = this->fadd.subspan(iP * Dimension::N, Dimension::N);
        auto grad     = this->grad.subspan(iP * Dimension::N, Dimension::N);
        auto ForceMat = this->ForceMat.subspan(iP * Dimension::NFF, Dimension::NFF);
        auto EMat     = this->EMat.subspan(iP * Dimension::FF, Dimension::FF);
        auto T        = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        auto V        = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        auto vpes     = this->vpes.subspan(iP, 1);
        auto alpha    = this->alpha.subspan(iP, 1);

        /////////////////////////////////////////////////////////////////
        // smooth dynamics force
        if (NAForce_type == NAForcePolicy::BOSD || NAForce_type == NAForcePolicy::NAFSD) {
            Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                             Kernel_Representation::inp_repr_type,    //
                                             Kernel_Representation::nuc_repr_type,    //
                                             SpacePolicy::L);
            ARRAY_CLEAR(fadd.data(), Dimension::N);  // additional force
            elec_utils::calc_distorted_rho(wrho.data(), rho_ele.data(), 1, 0, alpha[0]);
            double Ew_old = elec_utils::calc_ElectricalEnergy(EMat.data(), wrho.data(), occ_nuc[0]);
            alpha[0]      = calc_alpha(V.data());
            elec_utils::calc_distorted_rho(wrho.data(), rho_ele.data(), 1, 0, alpha[0]);
            double Ew_new = elec_utils::calc_ElectricalEnergy(EMat.data(), wrho.data(), occ_nuc[0]);
            elec_utils::calc_distorted_force(ftmp.data(), EMat.data(), ForceMat.data(), wrho.data(), rho_ele.data(),
                                             alpha[0]);

            // non-linear force
            for (int j = 0; j < Dimension::N; ++j) fadd[j] += ftmp[j];
            Epot[0] = vpes[0] + Ew_new;

            // region force
            if (Ew_new != Ew_old) {
                double vdotd = 0.0e0;
                for (int j = 0; j < Dimension::N; ++j) vdotd += p[j] / m[j] * ForceMat[j * Dimension::FF + 1];
                double xsolve = (Ew_new - Ew_old) / (dt_ptr[0]) / vdotd;  /// ??? scale
                for (int j = 0; j < Dimension::N; ++j) fadd[j] += xsolve * ForceMat[j * Dimension::FF + 1];
                for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = wrho[ik];
            }

            Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                             Kernel_Representation::nuc_repr_type,    //
                                             Kernel_Representation::inp_repr_type,    //
                                             SpacePolicy::L);
            Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                             Kernel_Representation::nuc_repr_type,    //
                                             Kernel_Representation::inp_repr_type,    //
                                             SpacePolicy::L);
        }

        switch (NAForce_type) {
            case NAForcePolicy::BO: {
                for (int j = 0, idxdV0 = 0; j < Dimension::N; ++j, idxdV0 += Dimension::FF)
                    f[j] = ForceMat[j * Dimension::FF + (occ_nuc[0]) * Dimension::Fadd1];
                break;
            }
            case NAForcePolicy::EHR: {
                Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {  // for both dV & dE (only for FMO-like model)
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        auto   Forceb0 = ForceMat.subspan(b0FF, Dimension::FF);
                        double fb0     = std::real(ARRAY_TRACE2(rho_nuc.data(), Forceb0.data(),  //
                                                            Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj] = fb0 * ForceMat[bjbb] / ForceMat[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        auto dVj = ForceMat.subspan(jFF, Dimension::FF);
                        f[j]     = std::real(ARRAY_TRACE2(rho_nuc.data(), dVj.data(),  //
                                                      Dimension::F, Dimension::F));
                    }
                }
                Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 SpacePolicy::L);
                break;
            }
            case NAForcePolicy::NAF: {
                Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 SpacePolicy::L);
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        auto   Forceb0 = ForceMat.subspan(b0FF, Dimension::FF);
                        double fb0     = Forceb0[(occ_nuc[0]) * Dimension::Fadd1];
                        double fprojb0 = std::real(ARRAY_TRACE2_OFFD(rho_nuc.data(), Forceb0.data(),  //
                                                                     Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj]     = fb0 * ForceMat[bjbb] / ForceMat[b0bb];  // @bug
                            fproj[bj] = fprojb0 * ForceMat[bjbb] / ForceMat[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        auto dVj = ForceMat.subspan(jFF, Dimension::FF);
                        f[j]     = dVj[(occ_nuc[0]) * Dimension::Fadd1];
                        fproj[j] = std::real(ARRAY_TRACE2_OFFD(rho_nuc.data(), dVj.data(), Dimension::F, Dimension::F));
                    }
                }
                if (offd_projected) {  // then the offdiagonal force is projected
                    double fdotv = 0.0e0, Pdotv = 0.0e0;
                    for (int j = 0; j < Dimension::N; ++j) fdotv += fproj[j] * p[j] / m[j], Pdotv += p[j] * p[j] / m[j];
                    if (Pdotv > 1.0E-8) {
                        for (int j = 0; j < Dimension::N; ++j) fproj[j] -= fdotv / Pdotv * p[j];
                    }
                }
                for (int j = 0; j < Dimension::N; ++j) f[j] += fproj[j];
                Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 SpacePolicy::L);
                break;
            }
            case NAForcePolicy::BOSD:
            case NAForcePolicy::NAFSD: {                                                  // smooth dynamics
                Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 SpacePolicy::L);

                // calculation of additional force
                ARRAY_CLEAR(fadd.data(), Dimension::N);
                elec_utils::calc_distorted_rho(wrho.data(), rho_ele.data(), 1, 0, alpha[0]);
                double Ew_old = elec_utils::calc_ElectricalEnergy(EMat.data(), wrho.data(), occ_nuc[0]);
                alpha[0]      = calc_alpha(V.data());
                elec_utils::calc_distorted_rho(wrho.data(), rho_ele.data(), 1, 0, alpha[0]);
                double Ew_new = elec_utils::calc_ElectricalEnergy(EMat.data(), wrho.data(), occ_nuc[0]);
                elec_utils::calc_distorted_force(ftmp.data(), EMat.data(), ForceMat.data(), wrho.data(), rho_ele.data(),
                                                 alpha[0]);

                // add non-linear force
                for (int j = 0; j < Dimension::N; ++j) fadd[j] += ftmp[j];
                Epot[0] = vpes[0] + Ew_new;

                // add regional force
                if (Ew_new != Ew_old) {
                    double vdotd = 0.0e0;
                    for (int j = 0; j < Dimension::N; ++j) vdotd += p[j] / m[j] * ForceMat[j * Dimension::FF + 1];
                    double xsolve = (Ew_new - Ew_old) / (dt_ptr[0]) / vdotd;  // scale???
                    for (int j = 0; j < Dimension::N; ++j) fadd[j] += xsolve * ForceMat[j * Dimension::FF + 1];
                }
                for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = wrho[ik];  // @nuc_repr_type

                // add weighted force
                if (FORCE_OPT::BATH_FORCE_BILINEAR) {
                    int& B   = FORCE_OPT::nbath;
                    int& J   = FORCE_OPT::Nb;
                    int  JFF = J * Dimension::FF;
                    for (int b = 0, bj = 0, b0FF = 0, b0bb = 0; b < B;
                         ++b, b0FF += JFF, b0bb += (JFF + Dimension::Fadd1)) {
                        auto   Forceb0 = ForceMat.subspan(b0FF, Dimension::FF);
                        double fb0 =
                            std::real(ARRAY_TRACE2_DIAG(rho_nuc.data(), Forceb0.data(), Dimension::F, Dimension::F));
                        double fprojb0 =
                            std::real(ARRAY_TRACE2_OFFD(rho_nuc.data(), Forceb0.data(), Dimension::F, Dimension::F));
                        for (int j = 0, bjbb = b0bb; j < J; ++j, ++bj, bjbb += Dimension::FF) {
                            f[bj]     = fb0 * ForceMat[bjbb] / ForceMat[b0bb];  // @bug
                            fproj[bj] = fprojb0 * ForceMat[bjbb] / ForceMat[b0bb];
                        }
                    }
                } else {
                    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
                        auto dVj = ForceMat.subspan(jFF, Dimension::FF);
                        f[j]     = std::real(ARRAY_TRACE2_DIAG(rho_nuc.data(), dVj.data(), Dimension::F, Dimension::F));
                        fproj[j] = std::real(ARRAY_TRACE2_OFFD(rho_nuc.data(), dVj.data(), Dimension::F, Dimension::F));
                    }
                }
                if (NAForce_type == NAForcePolicy::NAFSD) {  // twice off-diagonal force? @bug? check it please
                    if (offd_projected) {
                        double fdotR = 0.0e0, PdotR = 0.0e0;
                        for (int j = 0; j < Dimension::N; ++j)
                            fdotR += fproj[j] * p[j] / m[j], PdotR += p[j] * p[j] / m[j];
                        if (PdotR > 1.0e-8) {
                            for (int j = 0; j < Dimension::N; ++j) fproj[j] += fproj[j] - fdotR / PdotR * p[j];
                        }
                    } else {
                        for (int j = 0; j < Dimension::N; ++j) fproj[j] += fproj[j];
                    }
                }
                for (int j = 0; j < Dimension::N; ++j) f[j] += fproj[j];
                Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 SpacePolicy::L);
                Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                                 Kernel_Representation::nuc_repr_type,    //
                                                 Kernel_Representation::inp_repr_type,    //
                                                 SpacePolicy::L);
                break;
            }
        }
        for (int j = 0; j < Dimension::N; ++j) f[j] += grad[j];  // electronic independent force
        for (int j = 0; j < Dimension::N; ++j) f[j] += fadd[j];  // additional force
    }
    return stat;
}

NAForcePolicy::_type Kernel_NAForce::NAForce_type = NAForcePolicy::EHR;

namespace FORCE_OPT {
int  nbath               = 1;
int  Nb                  = 1;
bool BATH_FORCE_BILINEAR = false;
};  // namespace FORCE_OPT

};  // namespace PROJECT_NS

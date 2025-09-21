#include "psnd/Kernel_ExactPropagator.h"

#include "psnd/Kernel_Elec_Utils.h"
#include "psnd/Kernel_Representation.h"
#include "psnd/hash_fnv1a.h"
#include "psnd/linalg.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_ExactPropagator::getName() { return "Kernel_ExactPropagator"; }

int Kernel_ExactPropagator::getType() const { return utils::hash("Kernel_ExactPropagator"); }

void Kernel_ExactPropagator::setInputParam_impl(std::shared_ptr<Param> PM) {}

void Kernel_ExactPropagator::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    f    = DS->def(DATA::integrator::f);
    fadd = DS->def(DATA::integrator::fadd);
    p    = DS->def(DATA::integrator::p);
    ve   = DS->def(DATA::integrator::ve);
    minv = DS->def(DATA::integrator::minv);

    mask    = DS->def(DATA::integrator::forceeval::mask);
    dmask   = DS->def(DATA::integrator::forceeval::dmask);
    T       = DS->def(DATA::model::rep::T);
    grad    = DS->def(DATA::model::grad);
    dV      = DS->def(DATA::model::dV);
    dE      = DS->def(DATA::model::rep::dE);
    nac     = DS->def(DATA::model::rep::nac);
    eig     = DS->def(DATA::model::rep::eig);
    rho_nuc = DS->def(DATA::integrator::rho_nuc);
    c       = DS->def(DATA::integrator::c);
    Ekin    = DS->def(DATA::integrator::Ekin);
    dt      = DS->def(DATA::control::dt);

    // 分配 B_vec - Dim: [N]
    B_vec = DS->def(VARIABLE<psnd_real>("integrator.B_vec", &Dimension::shape_N, "@"));

    // 分配其他必需的变量
    alpha = DS->def(DATA::integrator::alpha);
    fproj = DS->def(DATA::integrator::tmp::fproj);
    ftmp  = DS->def(DATA::integrator::tmp::ftmp);
    wrho  = DS->def(DATA::integrator::tmp::wrho);

    // 分配缺少的变量
    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    rho_ele = DS->def(DATA::integrator::rho_ele);
    m       = DS->def(DATA::integrator::m);
    V       = DS->def(DATA::model::V);
    vpes    = DS->def(DATA::model::vpes);
    Epot    = DS->def(DATA::integrator::Epot);

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

Status& Kernel_ExactPropagator::initializeKernel_impl(Status& stat) { return stat; }

Status& Kernel_ExactPropagator::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto occ_nuc  = this->occ_nuc.subspan(iP, 1);
        auto rho_ele  = this->rho_ele.subspan(iP * Dimension::FF, Dimension::FF);
        auto rho_nuc  = this->rho_nuc.subspan(iP * Dimension::FF, Dimension::FF);
        auto f        = this->f.subspan(iP * Dimension::N, Dimension::N);
        auto p        = this->p.subspan(iP * Dimension::N, Dimension::N);
        auto m        = this->m.subspan(iP * Dimension::N, Dimension::N);
        auto fadd     = this->fadd.subspan(iP * Dimension::N, Dimension::N);
        auto fproj    = this->fproj.subspan(iP * Dimension::N, Dimension::N);
        auto grad     = this->grad.subspan(iP * Dimension::N, Dimension::N);
        auto ForceMat = this->ForceMat.subspan(iP * Dimension::NFF, Dimension::NFF);
        auto EMat     = this->EMat.subspan(iP * Dimension::FF, Dimension::FF);
        auto T        = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        auto V        = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        auto vpes     = this->vpes.subspan(iP, 1);
        auto alpha    = this->alpha.subspan(iP, 1);

        Kernel_Representation::transform(rho_nuc.data(), T.data(), Dimension::F,  //
                                         Kernel_Representation::inp_repr_type,    //
                                         Kernel_Representation::nuc_repr_type,    //
                                         SpacePolicy::L);

        for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
            auto dVj = ForceMat.subspan(jFF, Dimension::FF);
            f[j]     = dVj[(occ_nuc[0]) * Dimension::Fadd1];
            fproj[j] = std::real(ARRAY_TRACE2_OFFD(rho_nuc.data(), dVj.data(), Dimension::F, Dimension::F));
        }
    }

    return stat;
}

}  // namespace PROJECT_NS

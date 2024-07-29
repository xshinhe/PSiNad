#include "kids/Kernel_Update_U.h"

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Monodromy.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_U::getName() { return "Kernel_Update_U"; }

int Kernel_Update_U::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_U::setInputParam_impl(std::shared_ptr<Param> PM) {
    only_adjust            = _param->get_bool({"solver.only_adjust"}, LOC(), false);
    enable_update_c        = _param->get_bool({"solver.enable_update_c"}, LOC(), true);
    enable_update_cset     = _param->get_bool({"solver.enable_update_cset"}, LOC(), true);
    enable_update_rho_ele  = _param->get_bool({"solver.enable_update_rho_ele"}, LOC(), true);
    enable_update_rho_nuc  = _param->get_bool({"solver.enable_update_rho_nuc"}, LOC(), true);
    enable_update_rho_dual = _param->get_bool({"solver.enable_update_rho_dual"}, LOC(), true);
}

void Kernel_Update_U::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    eig           = DS->def(DATA::model::rep::eig);
    T             = DS->def(DATA::model::rep::T);
    dE            = DS->def(DATA::model::rep::dE);
    lam           = DS->def(DATA::model::rep::lam);
    R             = DS->def(DATA::model::rep::R);
    U             = DS->def(DATA::integrator::U);
    Udt           = DS->def(DATA::integrator::Udt);
    mono          = DS->def(DATA::integrator::monodromy::mono);
    monodt        = DS->def(DATA::integrator::monodromy::monodt);
    MFFtmp1       = DS->def(DATA::integrator::monodromy::MFFtmp1);
    MFFtmp2       = DS->def(DATA::integrator::monodromy::MFFtmp2);
    invexpidiagdt = DS->def(DATA::integrator::tmp::invexpidiagdt);
    c             = DS->def(DATA::integrator::c);
    cset          = DS->def(DATA::integrator::cset);
    rho_ele       = DS->def(DATA::integrator::rho_ele);
    rho_nuc       = DS->def(DATA::integrator::rho_nuc);
    rho_dual      = DS->def(DATA::integrator::rho_dual);
    c_init        = DS->def(DATA::init::c);
    cset_init     = DS->def(DATA::init::cset);
    rho_ele_init  = DS->def(DATA::init::rho_ele);
    rho_nuc_init  = DS->def(DATA::init::rho_nuc);
    rho_dual_init = DS->def(DATA::init::rho_dual);
    T_init        = DS->def(DATA::init::T);
    dt            = DS->def(DATA::flowcontrol::dt);
}

Status& Kernel_Update_U::initializeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_complex* U = this->U + iP * Dimension::FF;
        ARRAY_EYE(U, Dimension::F);
    }
    return stat;
}

Status& Kernel_Update_U::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        // local variables for iP-th of swarm
        kids_real*    eig           = this->eig + iP * Dimension::F;
        kids_real*    T             = this->T + iP * Dimension::FF;
        kids_real*    lam           = this->lam + iP * Dimension::F;
        kids_complex* R             = this->R + iP * Dimension::FF;
        kids_complex* U             = this->U + iP * Dimension::FF;
        kids_complex* Udt           = this->Udt + iP * Dimension::FF;
        kids_complex* c             = this->c + iP * Dimension::F;
        kids_complex* c_init        = this->c_init + iP * Dimension::F;
        kids_complex* cset          = this->cset + iP * Dimension::F;
        kids_complex* cset_init     = this->cset_init + iP * Dimension::F;
        kids_complex* rho_ele       = this->rho_ele + iP * Dimension::FF;
        kids_complex* rho_ele_init  = this->rho_ele_init + iP * Dimension::FF;
        kids_complex* rho_nuc       = this->rho_nuc + iP * Dimension::FF;
        kids_complex* rho_nuc_init  = this->rho_nuc_init + iP * Dimension::FF;
        kids_complex* rho_dual      = this->rho_dual + iP * Dimension::FF;
        kids_complex* rho_dual_init = this->rho_dual_init + iP * Dimension::FF;
        kids_real*    T_init        = this->T_init + iP * Dimension::FF;

        /**
         * Update propagator U
         */
        switch (Kernel_Representation::ele_repr_type) {
            case RepresentationPolicy::Diabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * eig[i] * scale * dt[0]);
                ARRAY_MATMUL3_TRANS2(Udt, T, invexpidiagdt, T, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            // DiabaticComplex: ...
            case RepresentationPolicy::Adiabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * lam[i] * scale * dt[0]);
                ARRAY_MATMUL3_TRANS2(Udt, R, invexpidiagdt, R, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            default:  // representation_policy::force, representation_policy::density
                throw kids_error("Unsupport Representation");
                break;
        }

        /**
         * Update c, cset, rho_ele, rho_nuc etc. (always keep synchronized with propagator U)
         */
        if (enable_update_c) {
            for (int i = 0; i < Dimension::F; ++i) c[i] = c_init[i];
            Kernel_Representation::transform(c, T_init, Dimension::F,               //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::ele_repr_type,  //
                                             SpacePolicy::H);
            ARRAY_MATMUL(c, U, c, Dimension::F, Dimension::F, 1);
            Kernel_Representation::transform(c, T, Dimension::F,                    //
                                             Kernel_Representation::ele_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::H);
        } else {
            enable_update_rho_ele = true;  // if not enable update of c; update of rho_ele must be enable
        }
        if (enable_update_rho_ele) {
            for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];
            Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::ele_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                             Kernel_Representation::ele_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
        } else {
            elec_utils::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
        }
        if (true || enable_update_rho_nuc) {  // @TODO
            for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = rho_nuc_init[ik];
            Kernel_Representation::transform(rho_nuc, T_init, Dimension::F,         //
                                             Kernel_Representation::inp_repr_type,  //
                                             Kernel_Representation::ele_repr_type,  //
                                             SpacePolicy::L);
            ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                             Kernel_Representation::ele_repr_type,  //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
            // where rho_nuc = rho_ele - Gamma, where Gamma is not evolutionary (not unitory!!!)
            if (only_adjust) {
                for (int i = 0; i < Dimension::FF; ++i) rho_nuc[i] = rho_ele[i] + (rho_nuc_init[i] - rho_ele_init[i]);
            }
        }
    }
    // trace on monodromy
    if (Kernel_Monodromy::enable) update_monodromy();
    return stat;
}

void Kernel_Update_U::update_monodromy() {
    int N0   = 0;
    int N1   = Dimension::N;
    int N2   = Dimension::N + Dimension::F;
    int N3   = 2 * Dimension::N + Dimension::F;
    int N4   = 2 * Dimension::N + 2 * Dimension::F;
    int N4N4 = N4 * N4;
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real*    eig    = this->eig + iP * Dimension::F;
        kids_real*    T      = this->T + iP * Dimension::FF;
        kids_real*    dE     = this->dE + iP * Dimension::NFF;
        kids_real*    mono   = this->mono + iP * N4N4;
        kids_real*    monodt = this->monodt + iP * N4N4;
        kids_complex* Udt    = this->Udt + iP * Dimension::FF;
        kids_complex* c      = this->c + iP * Dimension::F;

        // N0-N1: x, N1-N2:x_ele; N2-N3:p; N3-N4:p_ele
        ARRAY_EYE(monodt, N4);
        for (int i = 0, ik = 0; i < Dimension::F; ++i) {
            for (int k = 0; k < Dimension::F; ++k, ++ik) {
                monodt[(N1 + i) * N4 + (N1 + k)] = std::real(Udt[ik]);
                monodt[(N1 + i) * N4 + (N3 + k)] = -std::imag(Udt[ik]);
                monodt[(N3 + i) * N4 + (N1 + k)] = std::imag(Udt[ik]);
                monodt[(N3 + i) * N4 + (N3 + k)] = std::real(Udt[ik]);
            }
        }
        kids_complex im = phys::math::im;
        for (int j = 0, jik = 0; j < Dimension::N; ++j) {
            kids_complex* dxjUdt = MFFtmp1;  // as workspace
            kids_complex* c_tmp  = MFFtmp2;  // as workspace
            for (int i = 0, ii = 0, ik = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1, ++ik, ++jik) {
                    dxjUdt[ik] = (i == k) ? -std::exp(-im * eig[i] * dt[0]) * im * dE[jik] * dt[0]
                                          : dE[jik] / (eig[k] - eig[i]) *
                                                (std::exp(-im * eig[k] * dt[0]) - std::exp(-im * eig[i] * dt[0]));
                }
            }
            ARRAY_MATMUL3_TRANS2(dxjUdt, T, dxjUdt, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            ARRAY_MATMUL(c_tmp, dxjUdt, c, Dimension::F, Dimension::F, 1);
            for (int i = 0; i < Dimension::F; ++i) {
                monodt[(N1 + i) * N4 + (N0 + j)] = std::real(c_tmp[i]);
                monodt[(N3 + i) * N4 + (N0 + j)] = std::imag(c_tmp[i]);
            }
        }
        ARRAY_MATMUL(mono, monodt, mono, N4, N4, N4);
    }
}

};  // namespace PROJECT_NS
#include "kids/Sampling_MC.h"

#include <algorithm>

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Sampling_MC::getName() { return "Sampling_MC"; }

int Sampling_MC::getType() const { return utils::hash(FUNCTION_NAME); }


void Sampling_MC::setInputParam_impl(std::shared_ptr<Param> PM) {
    dt            = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d);  //
    alpha0        = _param->get_real({"solver.alpha0"}, LOC(), 1.0f);                  //
    width_scaling = _param->get_real({"solver.width_scaling"}, LOC(), 1.0f);           //
    break_thres   = _param->get_real({"solver.break_thres"}, LOC(), 1.0f);             //
    P_used0       = _param->get_int({"solver.P_initial"}, LOC(), 1);                   //
    max_clone     = _param->get_int({"solver.max_clone"}, LOC(), 0);                   //
    gamma         = _param->get_real({"solver.gamma"}, LOC(), 0.0f);                   //
    xi            = 1 + Dimension::F * gamma;

    impl_type          = _param->get_int({"solver.impl_type"}, LOC(), 0);           //
    samp_type          = _param->get_int({"solver.samp_type"}, LOC(), 0);           //
    aset_type          = _param->get_int({"solver.aset_type"}, LOC(), 0);           //
    time_displace_step = _param->get_int({"solver.time_displace_step"}, LOC(), 1);  //
}

void Sampling_MC::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    alpha    = DS->def(DATA::integrator::alpha);  //? N size check again
    Xcoeff   = DS->def(DATA::integrator::Xcoeff);
    Acoeff   = DS->def(DATA::integrator::Acoeff);
    dtAcoeff = DS->def(DATA::integrator::dtAcoeff);
    Hcoeff   = DS->def(DATA::integrator::Hcoeff);
    Hbasis   = DS->def(DATA::integrator::Hbasis);
    UXdt     = DS->def(DATA::integrator::UXdt);
    UYdt     = DS->def(DATA::integrator::UYdt);
    rhored   = DS->def(DATA::integrator::rhored);
    rhored2  = DS->def(DATA::integrator::rhored2);
    rhored3  = DS->def(DATA::integrator::rhored3);

    Snuc     = DS->def(DATA::integrator::Snuc);
    Sele     = DS->def(DATA::integrator::Sele);
    S        = DS->def(DATA::integrator::S);
    invS     = DS->def(DATA::integrator::invS);
    dtlnSnuc = DS->def(DATA::integrator::dtlnSnuc);
    dtSele   = DS->def(DATA::integrator::dtSele);
    L        = DS->def(DATA::integrator::GWP::L);
    L1       = DS->def(DATA::integrator::GWP::L1);
    L2       = DS->def(DATA::integrator::GWP::L2);
    R        = DS->def(DATA::integrator::GWP::R);
    R1       = DS->def(DATA::integrator::GWP::R1);
    R2       = DS->def(DATA::integrator::GWP::R2);
    S1       = DS->def(DATA::integrator::GWP::S1);
    S1h      = DS->def(DATA::integrator::GWP::S1h);
    invS1h   = DS->def(DATA::integrator::GWP::invS1h);
    S2       = DS->def(DATA::integrator::GWP::S2);
    S2h      = DS->def(DATA::integrator::GWP::S2h);
    invS2h   = DS->def(DATA::integrator::GWP::invS2h);
    Sx       = DS->def(DATA::integrator::GWP::Sx);

    Ekin          = DS->def(DATA::integrator::Ekin);
    g             = DS->def(DATA::integrator::g);
    clone_account = DS->def(DATA::integrator::clone_account);
    // pf_cross      = DS->def<kids_bool>(DATA::integrator::pf_cross);

    //
    Udt     = DS->def(DATA::integrator::Udt);
    H       = DS->def(DATA::model::rep::H);
    vpes    = DS->def(DATA::model::vpes);
    grad    = DS->def(DATA::model::grad);
    V       = DS->def(DATA::model::V);
    dV      = DS->def(DATA::model::dV);
    eig     = DS->def(DATA::model::rep::eig);  // check change to eig
    T       = DS->def(DATA::model::rep::T);
    dE      = DS->def(DATA::model::rep::dE);
    x       = DS->def(DATA::integrator::x);
    p       = DS->def(DATA::integrator::p);
    m       = DS->def(DATA::integrator::m);
    f       = DS->def(DATA::integrator::f);
    c       = DS->def(DATA::integrator::c);
    rho_ele = DS->def(DATA::integrator::rho_ele);
    rho_nuc = DS->def(DATA::integrator::rho_nuc);

    w       = DS->def(DATA::integrator::w);
    U       = DS->def(DATA::integrator::U);
    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    c_init  = DS->def(DATA::init::c);
    T_init  = DS->def(DATA::init::T);

    fun_diag_F = DS->def(DATA::integrator::tmp::fun_diag_F);
    fun_diag_P = DS->def(DATA::integrator::tmp::fun_diag_P);
    MatR_PP    = DS->def(DATA::integrator::tmp::MatR_PP);
    MatC_PP    = DS->def(DATA::integrator::tmp::MatC_PP);
    I_PP       = DS->def(DATA::integrator::tmp::I_PP);
    Ubranch    = DS->def(DATA::integrator::Ubranch);

    x_last    = DS->def(DATA::last::x);
    p_last    = DS->def(DATA::last::p);
    grad_last = DS->def(DATA::last::grad);
    dV_last   = DS->def(DATA::last::dV);
    g_last    = DS->def(DATA::last::g);
    c_last    = DS->def(DATA::last::c);

    P_used_ptr = DS->def(DATA::integrator::P_used);
    norm_ptr   = DS->def(DATA::integrator::norm);
    veF        = DS->def(DATA::integrator::veF);
    ve         = DS->def(DATA::integrator::ve);
}

Status& Sampling_MC::executeKernel_impl(Status& stat) {
    // @begin debug
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto x = this->x.subspan(iP * Dimension::N, Dimension::N);
        auto p = this->p.subspan(iP * Dimension::N, Dimension::N);
        //     for (int j = 0; j < Dimension::N; ++j) {
        //         x[j] = iP * 0.02 * (iP % 2 - 0.5) + 0.1 * j;
        //         p[j] = -iP * 0.02 * (iP % 2 - 0.5) + 0.2 * j;
        //     }
    }
    // @end debug

    // assuming Nucl & Elec are already performed. the following will give revision

    if (samp_type == 1) {  // overlap re-sampling
        for (int iP = 0; iP < Dimension::P; ++iP) {
            auto x = this->x.subspan(iP * Dimension::N, Dimension::N);
            auto p = this->p.subspan(iP * Dimension::N, Dimension::N);
            auto c = this->c.subspan(iP * Dimension::N, Dimension::N);
            for (int j = 0; j < Dimension::N; ++j) {
                x[j] = this->x[j];
                p[j] = this->p[j];
            }
        }
    }
    if (samp_type == 2) {  // neighbourhood re-sampling
        for (int iP = 0; iP < Dimension::P; ++iP) {
            auto x = this->x.subspan(iP * Dimension::N, Dimension::N);
            auto p = this->p.subspan(iP * Dimension::N, Dimension::N);
            auto c = this->c.subspan(iP * Dimension::N, Dimension::N);
            if (iP > 0) {  // fluctuation
                for (int j = 0; j < Dimension::N; ++j) {
                    double randu;
                    Kernel_Random::rand_gaussian(&randu);
                    randu = sqrt(iP * iP % (j + 1));
                    x[j]  = this->x[j] + width_scaling * randu / sqrt(Dimension::N * alpha0);
                    Kernel_Random::rand_gaussian(&randu);
                    randu = sqrt(iP % (j + 1));
                    p[j]  = this->p[j] + width_scaling * randu * sqrt(alpha0 / Dimension::N);
                }
            }
        }
    }
    if (samp_type == 3) {  /// time displaced re-sampling
        _kmodel->executeKernel(stat);
        _krepr->executeKernel(stat);
        _kforce->executeKernel(stat);
        for (int iP = 1; iP < Dimension::P; ++iP) {
            auto x_now       = this->x.subspan(iP * Dimension::N, Dimension::N);
            auto p_now       = this->p.subspan(iP * Dimension::N, Dimension::N);
            auto f_now       = this->f.subspan(iP * Dimension::N, Dimension::N);
            auto U_now       = this->U.subspan(iP * Dimension::FF, Dimension::FF);
            auto c_now       = this->c.subspan(iP * Dimension::F, Dimension::F);
            auto rho_nuc_now = this->rho_nuc.subspan(iP * Dimension::FF, Dimension::FF);
            auto x_prev      = this->x.subspan(std::max({iP - 2, 0}) * Dimension::N, Dimension::N);
            auto p_prev      = this->p.subspan(std::max({iP - 2, 0}) * Dimension::N, Dimension::N);
            auto f_prev      = this->f.subspan(std::max({iP - 2, 0}) * Dimension::N, Dimension::N);
            auto U_prev      = this->U.subspan(std::max({iP - 2, 0}) * Dimension::FF, Dimension::FF);
            auto eig_now     = this->eig.subspan(iP * Dimension::F, Dimension::F);
            auto T_now       = this->T.subspan(iP * Dimension::FF, Dimension::FF);
            auto Udt_now     = this->Udt.subspan(iP * Dimension::FF, Dimension::FF);

            kids_real signdt = (iP % 2 == 0) ? dt : -dt;

            for (int j = 0; j < Dimension::N; ++j) {  //
                x_now[j] = x_prev[j], p_now[j] = p_prev[j], f_now[j] = f_prev[j];
            }
            for (int istep_displace = 0; istep_displace < time_displace_step; ++istep_displace) {
                for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
                for (int j = 0; j < Dimension::N; ++j) x_now[j] += p_now[j] / m[j] * signdt;
                _kmodel->executeKernel(stat);
                _krepr->executeKernel(stat);
                switch (Kernel_Representation::ele_repr_type) {
                    case RepresentationPolicy::Diabatic: {
                        for (int i = 0; i < Dimension::F; ++i) {  //
                            fun_diag_F[i] = exp(-phys::math::im * eig_now[i] * signdt);
                        }
                        ARRAY_MATMUL3_TRANS2(Udt_now.data(), T_now.data(), fun_diag_F.data(), T_now.data(),  //
                                             Dimension::F, Dimension::F, 0, Dimension::F);
                        break;
                    }
                }
                ARRAY_MATMUL(U_now.data(), Udt_now.data(), U_prev.data(), Dimension::F, Dimension::F, Dimension::F);
                ARRAY_MATMUL(c_now.data(), U_now.data(), c.data(), Dimension::F, Dimension::F, 1);
                elec_utils::ker_from_c(rho_nuc_now.data(), c_now.data(), xi, gamma, Dimension::F);
                _kforce->executeKernel(stat);
                for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
            }
            // PRINT_ARRAY(x_now, 1, Dimension::N);
            // PRINT_ARRAY(p_now, 1, Dimension::N);
        }

        for (int iP = 0; iP < Dimension::P; ++iP) {
            auto w       = this->w.subspan(iP, 1);
            auto U       = this->U.subspan(iP * Dimension::FF, Dimension::FF);
            auto occ_nuc = this->occ_nuc.subspan(iP, 1);
            w[0]         = 1.0e0;               ///< initial measure
            occ_nuc[0]   = occ0;                ///< initial occupation
            ARRAY_EYE(U.data(), Dimension::F);  ///< initial propagator reset to identity
        }
    }
    _kmodel->executeKernel(stat);
    _krepr->executeKernel(stat);

    for (int i = 0; i < Dimension::N; ++i) alpha[i] = alpha0;
    for (int a = 0; a < Dimension::P; ++a) g[a] = 0.0e0;
    if (aset_type == 0) {
        for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = (a == 0) ? 1.0e0 : 0.0e0;
    } else if (aset_type == 1) {
        for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = 1.0e0;  // @only for test
    }

    // normalization of A
    P_used        = P_used0;
    P_used_ptr[0] = P_used;
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;

    _dataset->def(DATA::init::x, x);
    _dataset->def(DATA::init::p, p);
    _dataset->def(DATA::init::c, c);
    _dataset->def(DATA::init::rho_ele, rho_ele);
    _dataset->def(DATA::init::rho_nuc, rho_nuc);
    _dataset->def(DATA::init::T, T);

    norm_ptr[0] = 1.0e0;
    // PRINT_ARRAY(Acoeff, 1, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) g_last[a] = g[a];
    for (int aj = 0; aj < Dimension::PN; ++aj) x_last[aj] = x[aj];
    for (int aj = 0; aj < Dimension::PN; ++aj) p_last[aj] = p[aj];
    for (int ai = 0; ai < Dimension::PF; ++ai) c_last[ai] = c[ai];
    return stat;
}

};  // namespace PROJECT_NS

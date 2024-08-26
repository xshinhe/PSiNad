#include "kids/Kernel_Elec_Functions.h"

#include <algorithm>

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

double phi(double lambda, double N0_max, int F) {
    double x     = lambda * N0_max / 2;
    double term1 = 1.0e0;
    {  // hypergeometric(double a, double b, double c, double x)
        const double TOLERANCE = 1.0e-10;
        int          a = 1, b = 1, c = 1 + F;
        double       t = a * b * x / c;
        term1 += t;
        int n = 1;
        while (abs(t) > TOLERANCE) {
            a++, b++, c++, n++;
            t *= a * b * x / c / n;
            term1 += t;
        }
    }
    int Fact = 1;
    for (int i = 1; i < F; ++i) Fact *= i;
    term1        = 2 * (F - 2) * (F - 1) * term1 / Fact / F;
    double term2 = (6 - 2 * F - 2 * x) / Fact;
    return pow(lambda, F) * pow(2 - 2 * x, 1 - F) * (term1 + term2);
}

const std::string Kernel_Elec_Functions::getName() { return "Kernel_Elec_Functions"; }

int Kernel_Elec_Functions::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Elec_Functions::setInputParam_impl(std::shared_ptr<Param> PM) {
    occ0 = _param->get_int({"model.occ", "solver.occ"}, LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Dimension::F) throw std::runtime_error("occ >= F");

    gamma1 = _param->get_real({"solver.gamma"}, LOC(), elec_utils::gamma_wigner(Dimension::F));
    if (gamma1 < -1.5) gamma1 = elec_utils::gamma_opt(Dimension::F);
    if (gamma1 < -0.5) gamma1 = elec_utils::gamma_wigner(Dimension::F);
    gamma2 = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1    = (1 + Dimension::F * gamma1);
    xi2    = (1 + Dimension::F * gamma2);

    use_cv = _param->get_bool({"solver.use_cv"}, LOC(), false);
}

void Kernel_Elec_Functions::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    rho_ele = DS->def(DATA::integrator::rho_ele);
    T       = DS->def(DATA::model::rep::T);
    occ_nuc = DS->def(DATA::integrator::occ_nuc);

    w    = DS->def(DATA::integrator::w);
    w_CC = DS->def(DATA::integrator::w_CC);
    w_CP = DS->def(DATA::integrator::w_CP);
    w_PP = DS->def(DATA::integrator::w_PP);
    w_AA = DS->def(DATA::integrator::w_AA);
    w_AD = DS->def(DATA::integrator::w_AD);
    w_DD = DS->def(DATA::integrator::w_DD);
    wz_A = DS->def(DATA::integrator::wz_A);
    wz_D = DS->def(DATA::integrator::wz_D);
    ww_A = DS->def(DATA::integrator::ww_A);
    ww_D = DS->def(DATA::integrator::ww_D);

    K0   = DS->def(DATA::integrator::K0);
    K1   = DS->def(DATA::integrator::K1);
    K2   = DS->def(DATA::integrator::K2);
    K1QA = DS->def(DATA::integrator::K1QA);
    K2QA = DS->def(DATA::integrator::K2QA);
    K1DA = DS->def(DATA::integrator::K1DA);
    K2DA = DS->def(DATA::integrator::K2DA);
    K1QD = DS->def(DATA::integrator::K1QD);
    K2QD = DS->def(DATA::integrator::K2QD);
    K1DD = DS->def(DATA::integrator::K1DD);
    K2DD = DS->def(DATA::integrator::K2DD);

    ////

    KSHA = DS->def(DATA::integrator::KSHA);
    KTWA = DS->def(DATA::integrator::KTWA);
    KTWD = DS->def(DATA::integrator::KTWD);

    ////

    sqcw   = DS->def(DATA::integrator::sqcw);
    trKTWA = DS->def(DATA::integrator::trKTWA);
    trKTWD = DS->def(DATA::integrator::trKTWD);

    OpA = DS->def(DATA::integrator::OpA);
    OpB = DS->def(DATA::integrator::OpB);

    (DS->def(DATA::parameter::gamma0))[0] = _param->get_real({"solver.gamma0"}, LOC(), 0.0e0);
    (DS->def(DATA::parameter::gamma1))[0] = gamma1;
    (DS->def(DATA::parameter::gamma2))[0] = gamma2;
    (DS->def(DATA::parameter::gamma3))[0] = _param->get_real({"solver.gamma3"}, LOC(), 1.0e0);
    (DS->def(DATA::parameter::gammaw))[0] = elec_utils::gamma_wigner(Dimension::F);
    (DS->def(DATA::parameter::gammar))[0] = elec_utils::gamma_opt(Dimension::F);
    (DS->def(DATA::parameter::xi0))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gamma0))[0]);
    (DS->def(DATA::parameter::xi1))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gamma1))[0]);
    (DS->def(DATA::parameter::xi2))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gamma2))[0]);
    (DS->def(DATA::parameter::xi3))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gamma3))[0]);
    (DS->def(DATA::parameter::xiw))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gammaw))[0]);
    (DS->def(DATA::parameter::xir))[0]    = (1 + Dimension::F * (DS->def(DATA::parameter::gammar))[0]);

    kids_real* Is = DS->def(DATA::parameter::Is);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* Is0 = Is + iP * Dimension::FF;
        for (int i = 0, ik = 0; i < Dimension::F; ++i)
            for (int k = 0; k < Dimension::F; ++k, ++ik) Is0[ik] = (i == k) ? 1.0e0 : 0.0e0;
    }
}

Status& Kernel_Elec_Functions::initializeKernel_impl(Status& stat) {
    ww_A_init = _dataset->def_complex("init.ww_A", ww_A, Dimension::P);  // @bug!!!
    ww_D_init = _dataset->def_complex("init.ww_D", ww_D, Dimension::P);
    T_init    = _dataset->def_real("init.T", T, Dimension::PFF);

    executeKernel_impl(stat);
    double unit = 1.0e0;
    _dataset->def_real("integrator.1", &unit);
    _dataset->def_real("init.1", &unit);
    _dataset->def_complex("init.K0", K0, Dimension::PFF);
    _dataset->def_complex("init.K1", K1, Dimension::PFF);
    _dataset->def_complex("init.K2", K2, Dimension::PFF);
    _dataset->def_complex("init.K1QA", K1QA, Dimension::PFF);
    _dataset->def_complex("init.K2QA", K2QA, Dimension::PFF);
    _dataset->def_complex("init.K1DA", K1DA, Dimension::PFF);
    _dataset->def_complex("init.K2DA", K2DA, Dimension::PFF);
    _dataset->def_complex("init.K1QD", K1QD, Dimension::PFF);
    _dataset->def_complex("init.K2QD", K2QD, Dimension::PFF);
    _dataset->def_complex("init.K1DD", K1DD, Dimension::PFF);
    _dataset->def_complex("init.K2DD", K2DD, Dimension::PFF);

    _dataset->def_complex("init.KSHA", KSHA, Dimension::PFF);
    _dataset->def_complex("init.KTWA", KTWA, Dimension::PFF);
    _dataset->def_complex("init.KTWD", KTWD, Dimension::PFF);

    _dataset->def_complex("init.w", w, Dimension::P);
    _dataset->def_complex("init.wz_A", wz_A, Dimension::P);
    _dataset->def_complex("init.wz_D", wz_D, Dimension::P);
    ww_A_init = _dataset->def_complex("init.ww_A", ww_A, Dimension::P);
    ww_D_init = _dataset->def_complex("init.ww_D", ww_D, Dimension::P);
    T_init    = _dataset->def_real("init.T", T, Dimension::PFF);
    return stat;
}

Status& Kernel_Elec_Functions::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_complex* wz_A    = this->wz_A + iP;
        kids_complex* wz_D    = this->wz_D + iP;
        kids_complex* ww_A    = this->ww_A + iP;
        kids_complex* ww_D    = this->ww_D + iP;
        kids_complex* rho_ele = this->rho_ele + iP * Dimension::FF;
        kids_complex* K0      = this->K0 + iP * Dimension::FF;
        kids_complex* K1      = this->K1 + iP * Dimension::FF;
        kids_complex* K2      = this->K2 + iP * Dimension::FF;
        kids_complex* K1QA    = this->K1QA + iP * Dimension::FF;
        kids_complex* K2QA    = this->K2QA + iP * Dimension::FF;
        kids_complex* K1DA    = this->K1DA + iP * Dimension::FF;
        kids_complex* K2DA    = this->K2DA + iP * Dimension::FF;
        kids_complex* K1QD    = this->K1QD + iP * Dimension::FF;
        kids_complex* K2QD    = this->K2QD + iP * Dimension::FF;
        kids_complex* K1DD    = this->K1DD + iP * Dimension::FF;
        kids_complex* K2DD    = this->K2DD + iP * Dimension::FF;
        kids_complex* KSHA    = this->KSHA + iP * Dimension::FF;
        kids_complex* KTWA    = this->KTWA + iP * Dimension::FF;
        kids_complex* KTWD    = this->KTWD + iP * Dimension::FF;

        kids_real* sqcw   = this->sqcw + iP * Dimension::F;  // losed identity of sqc
        kids_real* trKTWA = this->trKTWA + iP;               // losed identity of sqc
        kids_real* trKTWD = this->trKTWD + iP;               // losed identity of sqc

        kids_real* T       = this->T + iP * Dimension::FF;
        int*       occ_nuc = this->occ_nuc + iP;

        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);

        // 1) Adiabatic representation
        /// parameters, windows(K), weights(w)
        wz_A[0] = 1.0e0;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            if (i == occ0) continue;
            wz_A[0] *= std::abs(rho_ele[occ0 * Dimension::Fadd1] - rho_ele[ii]);
        }

        int    max_pop = elec_utils::max_choose(rho_ele);
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);

        // K1Q simplex quantization
        int act = ((use_fall) ? occ_nuc[0] : max_pop);
        elec_utils::ker_from_rho(K1QA, rho_ele, 1, 0, Dimension::F, true, act);
        ARRAY_MAT_DIAG(K1DA, K1QA, Dimension::F);

        ww_A[0] = 4.0 - 1.0 / (max_val * max_val);
        ww_A[0] = std::min({std::abs(ww_A[0]), std::abs(ww_A_init[0])});

        // K2Q cutoff quantization (w2-window)
        elec_utils::ker_from_rho(K2QA, rho_ele, 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QA[i * Dimension::Fadd1]  //
                = (std::abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        ARRAY_MAT_DIAG(K2DA, K2QA, Dimension::F);

        // kernel for FSSH
        elec_utils::ker_from_rho(KSHA, rho_ele, 1, 0, Dimension::F, true, occ_nuc[0]);

        // kernel for TW
        elec_utils::ker_binning(KTWA, rho_ele, SQCPolicy::TRI);
        trKTWA[0] = std::real(ARRAY_TRACE1(KTWA, Dimension::F, Dimension::F));

        Kernel_Representation::transform(K1QA, T, Dimension::F,                 //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2QA, T, Dimension::F,                 //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1DA, T, Dimension::F,                 //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2DA, T, Dimension::F,                 //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(KSHA, T, Dimension::F,                 //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        // 2) Diabatic representation
        Kernel_Representation::transform(rho_ele, T, Dimension::F,         //
                                         RepresentationPolicy::Adiabatic,  //
                                         RepresentationPolicy::Diabatic,   //
                                         SpacePolicy::L);

        elec_utils::ker_from_rho(K0, rho_ele, 1.0e0, 0.0e0, Dimension::F);
        elec_utils::ker_from_rho(K1, rho_ele, xi1, gamma1, Dimension::F);
        elec_utils::ker_from_rho(K2, rho_ele, xi2, gamma2, Dimension::F);

        wz_D[0] = 1.0e0;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            if (i == occ0) continue;
            wz_D[0] *= std::abs(rho_ele[occ0 * Dimension::Fadd1] - rho_ele[ii]);
        }

        max_pop = elec_utils::max_choose(rho_ele);
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);

        elec_utils::ker_from_rho(K1QD, rho_ele, 1, 0, Dimension::F, true, max_pop);
        ARRAY_MAT_DIAG(K1DD, K1QD, Dimension::F);

        elec_utils::ker_from_rho(K2QD, rho_ele, 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QD[i * Dimension::Fadd1] = (std::abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        // if (use_strange_win)  // something else
        //     elec_utils::calc_distorted_rho(K2QD, rho_ele, 1, 0, 0.2);
        ARRAY_MAT_DIAG(K2DD, K2QD, Dimension::F);

        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        if (check_cxs) {
            double y = max_val - 0.5e0;
            ww_D[0]  = (23.0e0 / 2 * y - 72.0e0 * y * y + 140.0e0 * y * y * y + 480.0e0 * y * y * y * y -
                       840.0e0 * y * y * y * y * y) +
                      (9.0e0 / 4.0 - 9.0e0 * y * y - 420.0e0 * y * y * y * y + 1680.0e0 * y * y * y * y * y * y) *
                          atan(2.0e0 * y);
        }
        ww_D[0] = std::min({std::abs(ww_D[0]), std::abs(ww_D_init[0])});

        // general squeezed sqc
        if (false) {
            Kernel_Representation::transform(rho_ele_init, T_init, Dimension::F,    //
                                             Kernel_Representation::inp_repr_type,  //
                                             RepresentationPolicy::Diabatic,        //
                                             SpacePolicy::L);
            double vmax0 = 0.0e0, vsec0 = 0.0e0;
            double vmaxt = 0.0e0, vsect = 0.0e0;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                if (std::abs(rho_ele_init[ii]) > vmax0) {
                    vsec0 = vmax0;
                    vmax0 = std::abs(rho_ele_init[ii]);
                } else if (std::abs(rho_ele_init[ii]) > vsec0) {
                    vsec0 = std::abs(rho_ele_init[ii]);
                }
                if (std::abs(rho_ele[ii]) > vmax0) {
                    vsect = vmaxt;
                    vmaxt = std::abs(rho_ele[ii]);
                } else if (std::abs(rho_ele[ii]) > vsect) {
                    vsect = std::abs(rho_ele[ii]);
                }
            }

            double lambda1 = std::min({1 / vsect, 2 / (vmax0 + vsec0)});
            double lambda2 = std::max({1 / vmaxt, 1 / vmax0});
            if (lambda1 > lambda2) {
                ww_D[0] = phi(lambda1, vmax0, Dimension::F) - phi(lambda2, vmax0, Dimension::F);
            } else {
                ww_D[0] = 0.0e0;
            }
            Kernel_Representation::transform(rho_ele_init, T_init, Dimension::F,    //
                                             RepresentationPolicy::Diabatic,        //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
        }

        elec_utils::ker_binning(KTWD, rho_ele, SQCPolicy::TRI);
        if (count_exec <= 0) {  // only count at the beginning
            for (int k = 0; k < Dimension::F; ++k) {
                double radius = std::abs(2.0e0 - rho_ele[occ0 * Dimension::Fadd1] - rho_ele[k * Dimension::Fadd1]);
                sqcw[k]       = pow(radius, 3 - Dimension::F);
            }
        }
        if (sqc_init == 2) {  // overload for KTWD by shangyouhao
            int    imax = elec_utils::max_choose(rho_ele);
            double vmax = std::abs(rho_ele[imax * Dimension::Fadd1]);
            for (int ik = 0; ik < Dimension::FF; ++ik) KTWD[ik] = 0.0e0;
            if (vmax * vmax * 8.0e0 / 7.0e0 * (Dimension::F + 0.5e0) > 1) KTWD[imax * Dimension::Fadd1] = 1.0e0;
        }
        trKTWD[0] = std::real(ARRAY_TRACE1(KTWD, Dimension::F, Dimension::F));

        // 5) transform back from tcf_repr => inp_repr
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K0, T, Dimension::F,                   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1, T, Dimension::F,                   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2, T, Dimension::F,                   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1QD, T, Dimension::F,                 //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2QD, T, Dimension::F,                 //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1DD, T, Dimension::F,                 //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2DD, T, Dimension::F,                 //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
    }
    return stat;
}

};  // namespace PROJECT_NS

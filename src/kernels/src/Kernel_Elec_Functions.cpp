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

    auto Is = DS->def(DATA::parameter::Is);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto Is0 = Is.subspan(iP * Dimension::FF, Dimension::FF);
        for (int i = 0, ik = 0; i < Dimension::F; ++i)
            for (int k = 0; k < Dimension::F; ++k, ++ik) Is0[ik] = (i == k) ? 1.0e0 : 0.0e0;
    }
}

Status& Kernel_Elec_Functions::initializeKernel_impl(Status& stat) {
    if (_param->get_bool({"restart"}, LOC(), false)) {  //
        std::string loadfile = _param->get_string({"load"}, LOC(), "NULL");
        if (loadfile == "NULL" || loadfile == "" || loadfile == "null") loadfile = "restart.ds";
        _dataset->def(VARIABLE<kids_real>("integrator.1", &Dimension::shape_1, "@"))[0] = 1;
        _dataset->def(VARIABLE<kids_real>("init.1", &Dimension::shape_1, "@"))[0]       = 1;
        _dataset->def(VARIABLE<kids_complex>("init.K0", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K1", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K2", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K1QA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K2QA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K1DA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K2DA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K1QD", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K2QD", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K1DD", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.K2DD", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.KSHA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.KTWA", &Dimension::shape_PFF, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.KTWD", &Dimension::shape_PFF, "@"), loadfile);

        _dataset->def(VARIABLE<kids_complex>("init.w", &Dimension::shape_P, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.wz_A", &Dimension::shape_P, "@"), loadfile);
        _dataset->def(VARIABLE<kids_complex>("init.wz_D", &Dimension::shape_P, "@"), loadfile);

        // ? check if it is used ?
        ww_A_init = _dataset->def(VARIABLE<kids_complex>("init.ww_A", &Dimension::shape_P, "@"), loadfile);
        ww_D_init = _dataset->def(VARIABLE<kids_complex>("init.ww_D", &Dimension::shape_P, "@"), loadfile);
        T_init    = _dataset->def(VARIABLE<kids_real>("init.T", &Dimension::shape_PFF, "@"), loadfile);
        return stat;
    }

    ww_A_init = span<kids_complex>();  // blank for asking initialization
    ww_D_init = span<kids_complex>();  // blank for asking initialization
    executeKernel_impl(stat);
    _dataset->def(VARIABLE<kids_real>("integrator.1", &Dimension::shape_1, "@"))[0] = 1;
    _dataset->def(VARIABLE<kids_real>("init.1", &Dimension::shape_1, "@"))[0]       = 1;
    _dataset->def(VARIABLE<kids_complex>("init.K0", &Dimension::shape_PFF, "@"), K0);
    _dataset->def(VARIABLE<kids_complex>("init.K1", &Dimension::shape_PFF, "@"), K1);
    _dataset->def(VARIABLE<kids_complex>("init.K2", &Dimension::shape_PFF, "@"), K2);
    _dataset->def(VARIABLE<kids_complex>("init.K1QA", &Dimension::shape_PFF, "@"), K1QA);
    _dataset->def(VARIABLE<kids_complex>("init.K2QA", &Dimension::shape_PFF, "@"), K2QA);
    _dataset->def(VARIABLE<kids_complex>("init.K1DA", &Dimension::shape_PFF, "@"), K1DA);
    _dataset->def(VARIABLE<kids_complex>("init.K2DA", &Dimension::shape_PFF, "@"), K2DA);
    _dataset->def(VARIABLE<kids_complex>("init.K1QD", &Dimension::shape_PFF, "@"), K1QD);
    _dataset->def(VARIABLE<kids_complex>("init.K2QD", &Dimension::shape_PFF, "@"), K2QD);
    _dataset->def(VARIABLE<kids_complex>("init.K1DD", &Dimension::shape_PFF, "@"), K1DD);
    _dataset->def(VARIABLE<kids_complex>("init.K2DD", &Dimension::shape_PFF, "@"), K2DD);
    _dataset->def(VARIABLE<kids_complex>("init.KSHA", &Dimension::shape_PFF, "@"), KSHA);
    _dataset->def(VARIABLE<kids_complex>("init.KTWA", &Dimension::shape_PFF, "@"), KTWA);
    _dataset->def(VARIABLE<kids_complex>("init.KTWD", &Dimension::shape_PFF, "@"), KTWD);
    _dataset->def(VARIABLE<kids_complex>("init.w", &Dimension::shape_P, "@"), w);
    _dataset->def(VARIABLE<kids_complex>("init.wz_A", &Dimension::shape_P, "@"), wz_A);
    _dataset->def(VARIABLE<kids_complex>("init.wz_D", &Dimension::shape_P, "@"), wz_D);
    // fetch initialization
    ww_A_init = _dataset->def(VARIABLE<kids_complex>("init.ww_A", &Dimension::shape_P, "@"), ww_A);
    ww_D_init = _dataset->def(VARIABLE<kids_complex>("init.ww_D", &Dimension::shape_P, "@"), ww_D);
    T_init    = _dataset->def(VARIABLE<kids_real>("init.T", &Dimension::shape_PFF, "@"), T);
    return stat;
}

Status& Kernel_Elec_Functions::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto wz_A    = this->wz_A.subspan(iP, 1);
        auto wz_D    = this->wz_D.subspan(iP, 1);
        auto ww_A    = this->ww_A.subspan(iP, 1);
        auto ww_D    = this->ww_D.subspan(iP, 1);
        auto rho_ele = this->rho_ele.subspan(iP * Dimension::FF, Dimension::FF);
        auto K0      = this->K0.subspan(iP * Dimension::FF, Dimension::FF);
        auto K1      = this->K1.subspan(iP * Dimension::FF, Dimension::FF);
        auto K2      = this->K2.subspan(iP * Dimension::FF, Dimension::FF);
        auto K1QA    = this->K1QA.subspan(iP * Dimension::FF, Dimension::FF);
        auto K2QA    = this->K2QA.subspan(iP * Dimension::FF, Dimension::FF);
        auto K1DA    = this->K1DA.subspan(iP * Dimension::FF, Dimension::FF);
        auto K2DA    = this->K2DA.subspan(iP * Dimension::FF, Dimension::FF);
        auto K1QD    = this->K1QD.subspan(iP * Dimension::FF, Dimension::FF);
        auto K2QD    = this->K2QD.subspan(iP * Dimension::FF, Dimension::FF);
        auto K1DD    = this->K1DD.subspan(iP * Dimension::FF, Dimension::FF);
        auto K2DD    = this->K2DD.subspan(iP * Dimension::FF, Dimension::FF);
        auto KSHA    = this->KSHA.subspan(iP * Dimension::FF, Dimension::FF);
        auto KTWA    = this->KTWA.subspan(iP * Dimension::FF, Dimension::FF);
        auto KTWD    = this->KTWD.subspan(iP * Dimension::FF, Dimension::FF);

        auto sqcw    = this->sqcw.subspan(iP * Dimension::F, Dimension::F);  // losed identity of sqc
        auto trKTWA  = this->trKTWA.subspan(iP, 1);                          // losed identity of sqc
        auto trKTWD  = this->trKTWD.subspan(iP, 1);                          // losed identity of sqc
        auto T       = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        auto occ_nuc = this->occ_nuc.subspan(iP, 1);


        Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                         Kernel_Representation::inp_repr_type,    //
                                         RepresentationPolicy::Adiabatic,         //
                                         SpacePolicy::L);

        // 1) Adiabatic representation
        /// parameters, windows(K), weights(w)
        wz_A[0] = 1.0e0;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            if (i == occ0) continue;
            wz_A[0] *= std::abs(rho_ele[occ0 * Dimension::Fadd1] - rho_ele[ii]);
        }

        int    max_pop = elec_utils::max_choose(rho_ele.data());
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);

        // K1Q simplex quantization
        int act = ((use_fall) ? occ_nuc[0] : max_pop);
        elec_utils::ker_from_rho(K1QA.data(), rho_ele.data(), 1, 0, Dimension::F, true, act);
        ARRAY_MAT_DIAG(K1DA.data(), K1QA.data(), Dimension::F);

        ww_A[0] = 4.0 - 1.0 / (max_val * max_val);
        if (ww_A_init.size() > 0) ww_A[0] = std::min({std::abs(ww_A[0]), std::abs(ww_A_init[0])});

        // K2Q cutoff quantization (w2-window)
        elec_utils::ker_from_rho(K2QA.data(), rho_ele.data(), 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QA[i * Dimension::Fadd1]  //
                = (std::abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        ARRAY_MAT_DIAG(K2DA.data(), K2QA.data(), Dimension::F);

        // kernel for FSSH
        elec_utils::ker_from_rho(KSHA.data(), rho_ele.data(), 1, 0, Dimension::F, true, occ_nuc[0]);

        // kernel for TW
        elec_utils::ker_binning(KTWA.data(), rho_ele.data(), ElectronicSamplingPolicy::SQCtri);
        trKTWA[0] = std::real(ARRAY_TRACE1(KTWA.data(), Dimension::F, Dimension::F));

        Kernel_Representation::transform(K1QA.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2QA.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1DA.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2DA.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(KSHA.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        // 2) Diabatic representation
        Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                         RepresentationPolicy::Adiabatic,         //
                                         RepresentationPolicy::Diabatic,          //
                                         SpacePolicy::L);

        elec_utils::ker_from_rho(K0.data(), rho_ele.data(), 1.0e0, 0.0e0, Dimension::F);
        elec_utils::ker_from_rho(K1.data(), rho_ele.data(), xi1, gamma1, Dimension::F);
        elec_utils::ker_from_rho(K2.data(), rho_ele.data(), xi2, gamma2, Dimension::F);

        wz_D[0] = 1.0e0;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            if (i == occ0) continue;
            wz_D[0] *= std::abs(rho_ele[occ0 * Dimension::Fadd1] - rho_ele[ii]);
        }

        max_pop = elec_utils::max_choose(rho_ele.data());
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);

        elec_utils::ker_from_rho(K1QD.data(), rho_ele.data(), 1, 0, Dimension::F, true, max_pop);
        ARRAY_MAT_DIAG(K1DD.data(), K1QD.data(), Dimension::F);

        elec_utils::ker_from_rho(K2QD.data(), rho_ele.data(), 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QD[i * Dimension::Fadd1] = (std::abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        // if (use_strange_win)  // something else
        //     elec_utils::calc_distorted_rho(K2QD, rho_ele, 1, 0, 0.2);
        ARRAY_MAT_DIAG(K2DD.data(), K2QD.data(), Dimension::F);

        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        if (check_cxs) {
            double y = max_val - 0.5e0;
            ww_D[0]  = (23.0e0 / 2 * y - 72.0e0 * y * y + 140.0e0 * y * y * y + 480.0e0 * y * y * y * y -
                       840.0e0 * y * y * y * y * y) +
                      (9.0e0 / 4.0 - 9.0e0 * y * y - 420.0e0 * y * y * y * y + 1680.0e0 * y * y * y * y * y * y) *
                          atan(2.0e0 * y);
        }
        if (ww_D_init.size() > 0) ww_D[0] = std::min({std::abs(ww_D[0]), std::abs(ww_D_init[0])});

        // general squeezed sqc
        if (false) {
            Kernel_Representation::transform(rho_ele_init.data(), T_init.data(), Dimension::F,  //
                                             Kernel_Representation::inp_repr_type,              //
                                             RepresentationPolicy::Diabatic,                    //
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
            Kernel_Representation::transform(rho_ele_init.data(), T_init.data(), Dimension::F,  //
                                             RepresentationPolicy::Diabatic,                    //
                                             Kernel_Representation::inp_repr_type,              //
                                             SpacePolicy::L);
        }

        elec_utils::ker_binning(KTWD.data(), rho_ele.data(), ElectronicSamplingPolicy::SQCtri);
        if (count_exec <= 0) {  // only count at the beginning
            for (int k = 0; k < Dimension::F; ++k) {
                double radius = std::abs(2.0e0 - rho_ele[occ0 * Dimension::Fadd1] - rho_ele[k * Dimension::Fadd1]);
                sqcw[k]       = pow(radius, 3 - Dimension::F);
            }
        }
        if (sqc_init == 2) {  // overload for KTWD by shangyouhao
            int    imax = elec_utils::max_choose(rho_ele.data());
            double vmax = std::abs(rho_ele[imax * Dimension::Fadd1]);
            for (int ik = 0; ik < Dimension::FF; ++ik) KTWD[ik] = 0.0e0;
            if (vmax * vmax * 8.0e0 / 7.0e0 * (Dimension::F + 0.5e0) > 1) KTWD[imax * Dimension::Fadd1] = 1.0e0;
        }
        trKTWD[0] = std::real(ARRAY_TRACE1(KTWD.data(), Dimension::F, Dimension::F));

        // 5) transform back from tcf_repr => inp_repr
        Kernel_Representation::transform(rho_ele.data(), T.data(), Dimension::F,  //
                                         RepresentationPolicy::Diabatic,          //
                                         Kernel_Representation::inp_repr_type,    //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K0.data(), T.data(), Dimension::F,     //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1.data(), T.data(), Dimension::F,     //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2.data(), T.data(), Dimension::F,     //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1QD.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2QD.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K1DD.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2DD.data(), T.data(), Dimension::F,   //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
    }
    return stat;
}

};  // namespace PROJECT_NS

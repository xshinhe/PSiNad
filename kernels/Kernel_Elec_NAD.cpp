#include "Kernel_Elec_NAD.h"

#include "../core/linalg.h"
#include "../core/vars_list.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec_Utils.h"
// #include "elec_utils.h"
#include "Kernel_NADForce.h"
#include "Kernel_Random.h"
#include "Kernel_Representation.h"

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

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

namespace PROJECT_NS {

double calc_alpha(kids_real* V, int i = 0, int k = 1, int F = 2) {  // acoording to mix angle
    int ii = i * (F + 1), kk = k * (F + 1), ik = i * F + k;
    if (V[ii] == V[kk]) return 1.0e0;
    double res = 2.0e0 / phys::math::pi * atan(2 * V[ik] / (V[ii] - V[kk]));
    if (abs(res) < 1.0e-8) res = copysign(1.0e-8, res);
    return res;
}

int calc_wrho(kids_complex* wrho,   // distorted rho
              kids_complex* rho,    // rho_ele
              double        xi,     // xi must be 1
              double        gamma,  // gamma must be 0
              double        alpha) {
    // initialize distorted-density
    kids_real    L    = 1.0e0 - log(abs(alpha));  // @NOTE
    kids_complex norm = 0.0e0;
    for (int i = 0, ik = 0; i < Dimension::F; ++i) {
        for (int k = 0; k < Dimension::F; ++k, ++ik) {
            if (i == k) {
                wrho[ik] = copysign(1.0, std::real(rho[ik])) * std::pow(std::abs(rho[ik]), L);
                norm += wrho[ik];
            } else {
                wrho[ik] = xi * rho[ik];
            }
        }
    }
    // normalization of distorted-density
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) wrho[ii] /= norm;
    return 0;
}

double calc_Ew(kids_real* E, kids_complex* wrho, int occ) {
    double Ecalc = 0.0e0;
    switch (Kernel_NADForce::NADForce_type) {
        case NADForcePolicy::BO:
        case NADForcePolicy::CV: {
            Ecalc = E[occ];
            break;
        }
        default: {  // EHR, MIX, SD (Eto == Efrom will skip hopping procedure)
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                Ecalc += std::real(wrho[ii] * E[i]);
            }
            break;
        }
    }
    return Ecalc;
}

/**
 * @brief      the force driven from the shape of distorted-density W(\rho)
 *
 * @param      f1     The result
 * @param      E      adiabatic PES
 * @param      dE     The gradients of adiabatic PES
 * @param      wrho   The distorted-density W(\rho)
 * @param      rho    The density rho
 * @param[in]  xi     affine coefficient
 * @param[in]  gamma  affine coefficient
 * @param[in]  alpha  The distortion parameter
 *
 * @return     { description_of_the_return_value }
 */
int calc_distforce(kids_real*    f1,    // to be calculated
                   kids_real*    E,     // (input)
                   kids_real*    dE,    // (input)
                   kids_complex* wrho,  // distorted rho
                   kids_complex* rho,   // rho_ele
                   double        alpha) {
    // initialize distorted-density
    kids_real L            = 1.0e0 - log(abs(alpha));  // @NOTE
    kids_real rate_default = (L == 1.0e0) ? 1.0e0 : 0.0e0;

    // averaged adiabatic energy by distorted-density
    double Ew = 0.0e0;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) Ew += std::real(wrho[ii] * E[i]);

    // distortion force (when L= 1 [EHR], the distortion force is zero)
    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
        kids_real* dEj = dE + jFF;
        f1[j]          = 0.0e0;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            double rate  = ((std::real(rho[ii]) == 0.0e0) ? rate_default : std::real(wrho[ii] / rho[ii]));
            double coeff = (E[i] - (E[i] - Ew) * L * rate);
            for (int k = 0; k < Dimension::F; ++k) {
                if (i == k) continue;
                f1[j] += coeff * std::real(dEj[i * Dimension::F + k] / (E[k] - E[i]) * wrho[k * Dimension::F + i] -
                                           dEj[k * Dimension::F + i] / (E[i] - E[k]) * wrho[i * Dimension::F + k]);
            }
        }
    }
    return 0;
}

int hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm,  //
                    kids_real Efrom, kids_real Eto, int from, int to, bool reflect) {
    if (to == from) return from;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    kids_real coeffa = 0.0f, coeffb = 0.0f, coeffc = Eto - Efrom;
    for (int i = 0; i < Dimension::N; ++i) {
        coeffa += 0.5f * direction[i] * direction[i] / nm[i];
        coeffb += np[i] / nm[i] * direction[i];
    }
    coeffb /= coeffa, coeffc /= coeffa;  // normalization for safety

    kids_real coeffd = coeffb * coeffb - 4 * coeffc;
    if (coeffd > 0) {
        kids_real x1 = 0.5f * (-coeffb + sqrt(coeffd)), x2 = 0.5f * (-coeffb - sqrt(coeffd));
        kids_real xx = (std::abs(x1) < std::abs(x2)) ? x1 : x2;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return to;
    } else if (reflect) {  // 2008Algorithm
        kids_real xx = -coeffb;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return from;
    } else {  // 1990Algorithm, do nothing
        return from;
    }
    return from;
}

void Kernel_Elec_NAD::read_param_impl(Param* PM) {
    cmsh_type = NADPolicy::_from(PM->get<std::string>("cmsh_flag", LOC(), "CVSH"));

    alpha0 = PM->get<kids_real>("alpha0", LOC(), 0.5);
    gamma1 = PM->get<kids_real>("gamma", LOC(), elec_utils::gamma_wigner(Dimension::F));

    if (gamma1 < -1.5) gamma1 = elec_utils::gamma_opt(Dimension::F);
    if (gamma1 < -0.5) gamma1 = elec_utils::gamma_wigner(Dimension::F);

    gamma2          = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1             = (1 + Dimension::F * gamma1);
    xi2             = (1 + Dimension::F * gamma2);
    use_wmm         = PM->get<bool>("use_wmm", LOC(), false);
    use_sqc         = PM->get<bool>("use_sqc", LOC(), false);
    sqc_init        = PM->get<int>("sqc_init", LOC(), 0);
    use_fssh        = PM->get<bool>("use_fssh", LOC(), false);
    use_strange_win = PM->get<bool>("use_strange_win", LOC(), false);
    use_focus       = PM->get<bool>("use_focus", LOC(), false);
    use_fall        = PM->get<bool>("use_fall", LOC(), false);
    use_gdtwa       = PM->get<bool>("use_gdtwa", LOC(), false);
    only_adjust     = PM->get<bool>("only_adjust", LOC(), false);
    use_sum         = PM->get<bool>("use_sum", LOC(), false);
    check_cxs       = PM->get<bool>("check_cxs", LOC(), false);
    dt              = PM->get<double>("dt", LOC(), phys::time_d);

    cread_from_ds        = PM->get<bool>("cread_from_ds", LOC(), false);
    disable_inner_switch = PM->get<bool>("disable_inner_switch", LOC(), false);

    hopping_type1 = 0;
    hopping_type2 = 0;
    reflect       = true;
    use_cv        = false;
    dynamic_alpha = false;
    switch (cmsh_type) {
        case NADPolicy::EHR:
            Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;
            break;
        case NADPolicy::BOSH:
            Kernel_NADForce::NADForce_type = NADForcePolicy::BO;
            hopping_type1                  = 0;  // max(rho_ele)
            hopping_type2                  = 0;  // MASH direction
            reflect                        = true;
            break;
        case NADPolicy::CVSH:
            Kernel_NADForce::NADForce_type = NADForcePolicy::CV;
            hopping_type1                  = 1;  // max(rho_nuc)
            hopping_type2                  = 2;  // P direction
            reflect                        = false;
            use_cv                         = true;
            break;
        case NADPolicy::BOSD:
            Kernel_NADForce::NADForce_type = NADForcePolicy::BOSD;
            dynamic_alpha                  = true;
            break;
        case NADPolicy::CVSD:
            Kernel_NADForce::NADForce_type = NADForcePolicy::CVSD;
            dynamic_alpha                  = true;
            break;
    }
    hopping_type1 = PM->get<int>("hopping_type1", LOC(), hopping_type1);
    hopping_type2 = PM->get<int>("hopping_type2", LOC(), hopping_type2);
    use_cv        = PM->get<bool>("use_cv", LOC(), use_cv);
    reflect       = PM->get<bool>("reflect", LOC(), reflect);
    dynamic_alpha = PM->get<bool>("dynamic_alpha", LOC(), dynamic_alpha);
}

void Kernel_Elec_NAD::init_data_impl(DataSet* DS) {
    alpha                       = DS->def(DATA::integrator::alpha);
    Epot                        = DS->def(DATA::integrator::Epot);
    p                           = DS->def(DATA::integrator::p);
    m                           = DS->def(DATA::integrator::m);
    fadd                        = DS->def(DATA::integrator::fadd);
    ftmp                        = DS->def(DATA::integrator::tmp::ftmp);
    wrho                        = DS->def(DATA::integrator::tmp::wrho);
    vpes                        = DS->def(DATA::model::vpes);
    V                           = DS->def(DATA::model::V);
    E                           = DS->def(DATA::model::rep::E);
    dE                          = DS->def(DATA::model::rep::dE);
    T                           = DS->def(DATA::model::rep::T);
    H                           = DS->def(DATA::model::rep::H);
    direction                   = DS->def(DATA::integrator::tmp::direction);
    sqcw                        = DS->def(DATA::integrator::sqcw);
    sqcIA                       = DS->def(DATA::integrator::sqcIA);
    sqcID                       = DS->def(DATA::integrator::sqcID);
    at_samplingstep_finally_ptr = DS->def(DATA::iter::at_samplingstep_finally);
}

void Kernel_Elec_NAD::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_complex* w       = Kernel_Elec::w + iP;
        kids_complex* wz_A    = Kernel_Elec::wz_A + iP;
        kids_complex* wz_D    = Kernel_Elec::wz_D + iP;
        kids_complex* ww_A    = Kernel_Elec::ww_A + iP;
        kids_complex* ww_D    = Kernel_Elec::ww_D + iP;
        kids_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        kids_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        kids_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        kids_real*    T       = Kernel_Elec::T + iP * Dimension::FF;
        int*          occ_nuc = Kernel_Elec::occ_nuc + iP;
        kids_real*    alpha   = this->alpha + iP;
        kids_real*    V       = this->V + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////
        alpha[0] = (dynamic_alpha) ? calc_alpha(V) : alpha0;

        w[0] = kids_complex(Dimension::F);  ///< initial measure
        int iocc;
        Kernel_Random::rand_catalog(&iocc, 1, true, 0, Dimension::F - 1);
        iocc     = ((use_sum) ? iocc : Kernel_Elec::occ0);
        *occ_nuc = iocc;
        if (use_focus) {
            elec_utils::c_focus(c, xi1, gamma1, iocc, Dimension::F);
        } else if (use_sqc) {
            switch (sqc_init) {
                case 0: {  // traditional SQC
                    elec_utils::c_window(c, iocc, SQCPolicy::TRI, Dimension::F);
                    break;
                }
                case 1: {                                   // simplex SQC
                    elec_utils::c_sphere(c, Dimension::F);  ///< initial c on standard sphere
                    for (int i = 0; i < Dimension::F; ++i) c[i] = abs(c[i] * c[i]);
                    c[iocc] += 1.0e0;
                    for (int i = 0; i < Dimension::F; ++i) {
                        kids_real randu;
                        Kernel_Random::rand_uniform(&randu);
                        randu *= phys::math::twopi;
                        c[i] = sqrt(c[i]);
                        c[i] *= (cos(randu) + phys::math::im * sin(randu));
                    }
                    break;
                }
                case 2: {  // suggested by YHShang
                    elec_utils::c_sphere(c, Dimension::F);
                    break;
                }
                case 3: {  // suggested by Liu
                    elec_utils::c_window(c, iocc, SQCPolicy::TRI, Dimension::F);
                    break;
                }
                case 4: {  // suggested by Liu
                    elec_utils::c_window(c, iocc, SQCPolicy::TRI, Dimension::F);
                    double norm = 0.0e0;
                    for (int i = 0; i < Dimension::F; ++i) norm += std::abs(c[i] * c[i]);
                    xi1    = norm;
                    gamma1 = (xi1 - 1.0e0) / Dimension::F;
                    norm   = sqrt(norm);
                    for (int i = 0; i < Dimension::F; ++i) c[i] /= norm;
                    break;
                }
            }
        } else {
            elec_utils::c_sphere(c, Dimension::F);  ///< initial c on standard sphere
            // for (int i = 0; i < Dimension::F; ++i) c[i] = 1.0e0 * i;
        }
        if (cread_from_ds) {
            std::string init_nuclinp = _Param->get<std::string>("init_nuclinp", LOC());
            std::string open_file    = init_nuclinp;
            if (!isFileExists(init_nuclinp)) open_file = utils::concat(init_nuclinp, stat, ".ds");

            std::string   stmp, eachline;
            std::ifstream ifs(open_file);
            while (getline(ifs, eachline)) {
                if (eachline.find("init.c") != eachline.npos) {
                    getline(ifs, eachline);
                    for (int i = 0; i < Dimension::N; ++i) ifs >> c[i];
                }
            }
        }

        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele

        if (use_gdtwa) {
            xi1    = 1.0e0;
            gamma1 = 0.0e0;
            use_cv = false;

            double randu    = 1.0e0;
            double gamma_ou = phys::math::sqrthalf;
            double gamma_uu = 0.0e0;
            for (int j = 0; j < Dimension::F; ++j) {
                if (j == iocc) {
                    rho_ele[j * Dimension::Fadd1] = 1.0e0;
                    continue;
                }
                Kernel_Random::rand_uniform(&randu);
                randu                            = phys::math::halfpi * (int(randu / 0.25f) + 0.5);
                rho_ele[iocc * Dimension::F + j] = cos(randu) + phys::math::im * sin(randu);
                rho_ele[j * Dimension::F + iocc] = std::conj(rho_ele[iocc * Dimension::F + j]);
            }
            for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                for (int j = 0; j < Dimension::F; ++j, ++ij) {
                    if (i == iocc || j == iocc) continue;
                    rho_ele[ij] = rho_ele[iocc * Dimension::F + j] / rho_ele[iocc * Dimension::F + i];
                }
            }
            for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                for (int j = 0; j < Dimension::F; ++j, ++ij) {
                    if (i == j) {
                        rho_ele[ij] = (i == iocc) ? phys::math::iu : phys::math::iz;
                    } else if (i == iocc || j == iocc) {
                        rho_ele[ij] *= gamma_ou;
                    } else {
                        rho_ele[ij] *= gamma_uu;
                    }
                }
            }
        }

        if (use_sqc) {
            if (sqc_init == 0 || sqc_init == 1) {
                Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, 1.0, gamma1, Dimension::F, use_cv,
                                          iocc);  ///< initial rho_nuc
            } else if (sqc_init <= 3) {           // by youhao shang / Liu scheme
                double norm = 0.0;
                for (int i = 0; i < Dimension::F; ++i) norm += abs(rho_ele[i * Dimension::Fadd1]);
                Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi1 / norm, gamma1, Dimension::F, use_cv,
                                          iocc);  ///< initial rho_nuc
            } else if (sqc_init <= 4) {
                double norm = 0.0;
                for (int i = 0; i < Dimension::F; ++i) norm += abs(rho_ele[i * Dimension::Fadd1]);
                Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, 1.0, (norm - 1.0) / Dimension::F, Dimension::F, use_cv,
                                          iocc);  ///< initial rho_nuc
            }
        } else {
            Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv,
                                      iocc);  ///< initial rho_nuc
        }
        ARRAY_EYE(U, Dimension::F);  ///< initial propagator

        // BO occupation in adiabatic representation
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);
        *occ_nuc = elec_utils::max_choose(rho_nuc);
        if (use_fssh) *occ_nuc = elec_utils::pop_choose(rho_nuc);
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        // weight factor in tcf_repr
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);
        wz_A[0]        = std::abs(rho_ele[0] - rho_ele[3]);
        int    max_pop = elec_utils::max_choose(rho_ele);
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_A[0]        = 4.0 - 1.0 / (max_val * max_val);
        Kernel_Representation::transform(rho_ele, T, Dimension::F,         //
                                         RepresentationPolicy::Adiabatic,  //
                                         RepresentationPolicy::Diabatic,   //
                                         SpacePolicy::L);
        wz_D[0] = std::abs(rho_ele[0] - rho_ele[3]);
        max_pop = elec_utils::max_choose(rho_ele);
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
    }

    _DataSet->def("init.wz_A", Kernel_Elec::wz_A, Dimension::P);
    _DataSet->def("init.wz_D", Kernel_Elec::wz_D, Dimension::P);
    Kernel_Elec::ww_A_init    = _DataSet->def("init.ww_A", Kernel_Elec::ww_A, Dimension::P);
    Kernel_Elec::ww_D_init    = _DataSet->def("init.ww_D", Kernel_Elec::ww_D, Dimension::P);
    Kernel_Elec::c_init       = _DataSet->def("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->def("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::rho_nuc_init = _DataSet->def("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    Kernel_Elec::T_init       = _DataSet->def("init.T", Kernel_Elec::T, Dimension::PFF);

    at_samplingstep_finally_ptr[0] = true;
    exec_kernel(stat);
}

int Kernel_Elec_NAD::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        int*          occ_nuc      = Kernel_Elec::occ_nuc + iP;
        kids_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        kids_complex* c            = Kernel_Elec::c + iP * Dimension::F;
        kids_complex* c_init       = Kernel_Elec::c_init + iP * Dimension::F;
        kids_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        kids_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        kids_complex* rho_nuc_init = Kernel_Elec::rho_nuc_init + iP * Dimension::FF;
        kids_real*    T            = Kernel_Elec::T + iP * Dimension::FF;
        kids_real*    T_init       = Kernel_Elec::T_init + iP * Dimension::FF;
        kids_complex* K0           = Kernel_Elec::K0 + iP * Dimension::FF;
        kids_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        kids_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;
        kids_complex* K1QA         = Kernel_Elec::K1QA + iP * Dimension::FF;
        kids_complex* K2QA         = Kernel_Elec::K2QA + iP * Dimension::FF;
        kids_complex* K1DA         = Kernel_Elec::K1DA + iP * Dimension::FF;
        kids_complex* K2DA         = Kernel_Elec::K2DA + iP * Dimension::FF;
        kids_complex* K1QD         = Kernel_Elec::K1QD + iP * Dimension::FF;
        kids_complex* K2QD         = Kernel_Elec::K2QD + iP * Dimension::FF;
        kids_complex* K1DD         = Kernel_Elec::K1DD + iP * Dimension::FF;
        kids_complex* K2DD         = Kernel_Elec::K2DD + iP * Dimension::FF;
        kids_complex* ww_A_init    = Kernel_Elec::ww_A_init + iP;
        kids_complex* ww_D_init    = Kernel_Elec::ww_D_init + iP;
        kids_complex* ww_A         = Kernel_Elec::ww_A + iP;
        kids_complex* ww_D         = Kernel_Elec::ww_D + iP;

        kids_real*    alpha = this->alpha + iP;
        kids_real*    Epot  = this->Epot + iP;
        kids_real*    vpes  = this->vpes + iP;
        kids_real*    V     = this->V + iP * Dimension::FF;
        kids_real*    E     = this->E + iP * Dimension::F;
        kids_real*    dE    = this->dE + iP * Dimension::NFF;
        kids_real*    p     = this->p + iP * Dimension::N;
        kids_real*    m     = this->m + iP * Dimension::N;
        kids_real*    fadd  = this->fadd + iP * Dimension::N;
        kids_complex* H     = this->H + iP * Dimension::FF;

        //////////////////////////////////////////////////////////////////////

        // * additional evolution can be appended here

        // 1) transform from inp_repr => ele_repr
        for (int i = 0; i < Dimension::F; ++i) c[i] = c_init[i];  // @debug
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = rho_nuc_init[ik];
        Kernel_Representation::transform(c, T_init, Dimension::F,               //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::H);
        Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(rho_nuc, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);

        // 2) propagte along ele_repr
        ARRAY_MATMUL(c, U, c, 1, Dimension::F, Dimension::F);
        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);

        // 3) hopping in adiabatic representation
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);

        if (!disable_inner_switch) {
            // std::cout << "****\n";
            // step 1: determine where to hop (BOSH & CVSH)
            /// 1.1 calc Efrom
            kids_real Efrom, Eto;
            Efrom = calc_Ew(E, rho_nuc, *occ_nuc);
            /// 1.2 calc to
            int to;
            switch (hopping_type1) {
                case 0: {
                    to = elec_utils::max_choose(rho_ele);
                    break;
                }
                case 1: {
                    to = elec_utils::max_choose(rho_nuc);
                    break;
                }
                case 2: {
                    to = elec_utils::pop_choose(rho_ele);
                    break;
                }
                case 3: {
                    to = elec_utils::hopping_choose(rho_ele, H, *occ_nuc, scale * dt);
                    break;
                }
                case 4: {
                    to = elec_utils::pop_choose(rho_nuc);
                    break;
                }
                case 5: {
                    to = elec_utils::pop_neg_choose(rho_nuc);
                    break;
                }
            }
            Eto = calc_Ew(E, rho_nuc, to);
            // step 2: determine direction to hop
            switch (hopping_type2) {
                case 0: {
                    elec_utils::hopping_direction(direction, E, dE, rho_ele, *occ_nuc, to);
                    break;
                }
                case 1: {
                    elec_utils::hopping_direction(direction, dE, *occ_nuc, to);
                    break;
                }
                case 2: {
                    for (int j = 0; j < Dimension::N; ++j) direction[j] = p[j];
                    break;
                }
                case 3: {
                    for (int j = 0; j < Dimension::N; ++j)
                        direction[j] = dE[j * Dimension::FF + to * Dimension::Fadd1] -
                                       dE[j * Dimension::FF + (*occ_nuc) * Dimension::Fadd1];
                    break;
                }
            }
            // step 3: try hop
            *occ_nuc = hopping_impulse(direction, p, m, Efrom, Eto, *occ_nuc, to, reflect);
            Epot[0]  = vpes[0] + ((*occ_nuc == to) ? Eto : Efrom);
        }

        // if (*occ_nuc == to && Eto != Efrom) std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!\n";

        // smooth dynamics (BOSD & CVSD)
        if (cmsh_type == NADPolicy::BOSD || cmsh_type == NADPolicy::CVSD) {
            ARRAY_CLEAR(fadd, Dimension::N);

            // Win is for calculate energy
            calc_wrho(wrho, rho_ele, xi1, gamma1, alpha[0]);
            double Ew_old = calc_Ew(E, wrho, *occ_nuc);
            calc_distforce(ftmp, E, dE, wrho, rho_ele, alpha[0]);

            if (dynamic_alpha) alpha[0] = calc_alpha(V);
            calc_wrho(wrho, rho_ele, xi1, gamma1, alpha[0]);
            double Ew_new = calc_Ew(E, wrho, *occ_nuc);
            calc_distforce(ftmp, E, dE, wrho, rho_ele, alpha[0]);
            // non-linear force
            for (int j = 0; j < Dimension::N; ++j) fadd[j] += ftmp[j];
            Epot[0] = vpes[0] + Ew_new;

            // region force
            if (Ew_new != Ew_old) {
                double vdotd = 0.0e0;
                for (int j = 0; j < Dimension::N; ++j) vdotd += p[j] / m[j] * dE[j * Dimension::FF + 1];
                double xsolve = (Ew_new - Ew_old) / (scale * dt) / vdotd;
                for (int j = 0; j < Dimension::N; ++j) fadd[j] += xsolve * dE[j * Dimension::FF + 1];
            }
            for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = wrho[ik];
        }

        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        if ((use_sqc && reflect) || (use_sqc && only_adjust)) {                     // only adjusted
            Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                             RepresentationPolicy::Adiabatic,       //
                                             Kernel_Representation::inp_repr_type,  //
                                             SpacePolicy::L);
            double norm = 0.0e0;
            for (int i = 0; i < Dimension::F; ++i) norm += abs(rho_ele[i * Dimension::Fadd1]);
            // Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi1/norm, gamma1, Dimension::F, use_cv,
            //                          iocc);  ///< initial rho_nuc /// adjust only support for sqc_init=0!!!
            for (int i = 0; i < Dimension::FF; ++i) rho_nuc[i] = rho_ele[i] + (rho_nuc_init[i] - rho_ele_init[i]);
            Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                             Kernel_Representation::inp_repr_type,  //
                                             RepresentationPolicy::Adiabatic,       //
                                             SpacePolicy::L);
        }

        if (!at_samplingstep_finally_ptr[0]) continue;

        // 4) calculated TCF in adiabatic rep & diabatic rep respectively
        // 4-1) Adiabatic rep
        int    max_pop = elec_utils::max_choose(rho_ele);
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_A[0]        = 4.0 - 1.0 / (max_val * max_val);
        ww_A[0]        = std::min({abs(ww_A[0]), abs(ww_A_init[0])});

        Kernel_Elec::ker_from_rho(K1QA, rho_ele, 1, 0, Dimension::F, true, ((use_fall) ? *occ_nuc : max_pop));
        Kernel_Elec::ker_from_rho(K2QA, rho_ele, 1, 0, Dimension::F, true, ((use_fall) ? *occ_nuc : max_pop));
        for (int i = 0; i < Dimension::F; ++i) {
            K2QA[i * Dimension::Fadd1] = (abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        if (use_fssh) { Kernel_Elec::ker_from_rho(K2QA, rho_ele, 1, 0, Dimension::F, true, *occ_nuc); }
        if (use_sqc) {
            elec_utils::ker_binning(K2QA, rho_ele, SQCPolicy::TRI);
            sqcIA[0] = 0;
            for (int i = 0; i < Dimension::F; ++i) sqcIA[0] += std::real(K2QA[i * Dimension::Fadd1]);
        }

        ARRAY_MAT_DIAG(K1DA, K1QA, Dimension::F);
        ARRAY_MAT_DIAG(K2DA, K2QA, Dimension::F);

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
        // 4-2) Diabatic rep
        Kernel_Representation::transform(rho_ele, T, Dimension::F,
                                         RepresentationPolicy::Adiabatic,  //
                                         RepresentationPolicy::Diabatic,   //
                                         SpacePolicy::L);

        Kernel_Elec::ker_from_rho(K0, rho_ele, 1, 0, Dimension::F);
        Kernel_Elec::ker_from_rho(K1, rho_ele, xi1, gamma1, Dimension::F);
        Kernel_Elec::ker_from_rho(K2, rho_ele, xi2, gamma2, Dimension::F);

        max_pop = elec_utils::max_choose(rho_ele);  // (in dia rep)
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);

        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        if (check_cxs) {
            double y = max_val - 0.5e0;
            ww_D[0]  = (23.0e0 / 2 * y - 72.0e0 * y * y + 140.0e0 * y * y * y + 480.0e0 * y * y * y * y -
                       840.0e0 * y * y * y * y * y) +
                      (9.0e0 / 4.0 - 9.0e0 * y * y - 420.0e0 * y * y * y * y + 1680.0e0 * y * y * y * y * y * y) *
                          atan(2.0e0 * y);
        }
        ww_D[0] = std::min({abs(ww_D[0]), abs(ww_D_init[0])});
        {                                                                           // general squeezed sqc
            Kernel_Representation::transform(rho_ele_init, T_init, Dimension::F,    //
                                             Kernel_Representation::inp_repr_type,  //
                                             RepresentationPolicy::Diabatic,        //
                                             SpacePolicy::L);
            double vmax0 = 0.0e0, vsec0 = 0.0e0;
            double vmaxt = 0.0e0, vsect = 0.0e0;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                if (abs(rho_ele_init[ii]) > vmax0) {
                    vsec0 = vmax0;
                    vmax0 = abs(rho_ele_init[ii]);
                } else if (abs(rho_ele_init[ii]) > vsec0) {
                    vsec0 = abs(rho_ele_init[ii]);
                }
                if (abs(rho_ele[ii]) > vmax0) {
                    vsect = vmaxt;
                    vmaxt = abs(rho_ele[ii]);
                } else if (abs(rho_ele[ii]) > vsect) {
                    vsect = abs(rho_ele[ii]);
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

        Kernel_Elec::ker_from_rho(K1QD, rho_ele, 1, 0, Dimension::F, true, max_pop);
        Kernel_Elec::ker_from_rho(K2QD, rho_ele, 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QD[i * Dimension::Fadd1] = (abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        if (use_strange_win) calc_wrho(K2QD, rho_ele, 1, 0, 0.2);
        if (use_sqc) {
            elec_utils::ker_binning(K2QD, rho_ele, SQCPolicy::TRI);
            if (count_exec <= 0) {  // only count at the beginning
                for (int k = 0; k < Dimension::F; ++k) {
                    double radius =
                        abs(2.0e0 - rho_ele[Kernel_Elec::occ0 * Dimension::Fadd1] - rho_ele[k * Dimension::Fadd1]);
                    // radius  = sqrt(radius * radius + 0.0025e0); // @bad soft radius
                    sqcw[k] = pow(radius, 3 - Dimension::F);
                }
            }

            sqcID[0] = 0;
            for (int i = 0; i < Dimension::F; ++i) sqcID[0] += std::real(K2QD[i * Dimension::Fadd1]);

            if (sqc_init == 2) {  // overload for K2QD
                int    imax = elec_utils::max_choose(rho_ele);
                double vmax = std::abs(rho_ele[imax * Dimension::Fadd1]);
                for (int ik = 0; ik < Dimension::FF; ++ik) K2QD[ik] = 0.0e0;
                if (vmax * vmax * 8.0e0 / 7.0e0 * (Dimension::F + 0.5e0) > 1) K2QD[imax * Dimension::Fadd1] = 1.0e0;
            }
        }

        ARRAY_MAT_DIAG(K1DD, K1QD, Dimension::F);
        ARRAY_MAT_DIAG(K2DD, K2QD, Dimension::F);

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

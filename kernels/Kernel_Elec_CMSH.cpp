#include "Kernel_Elec_CMSH.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec_CMM.h"
#include "Kernel_Elec_MMSH.h"
#include "Kernel_Elec_SH.h"
#include "Kernel_Elec_SQC.h"
#include "Kernel_NADForce.h"
#include "Kernel_Random.h"
#include "Kernel_Representation.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

namespace PROJECT_NS {

double calc_alpha(num_real* V, int i = 0, int k = 1, int F = 2) {  // acoording to mix angle
    int ii = i * (F + 1), kk = k * (F + 1), ik = i * F + k;
    if (V[ii] == V[kk]) return 1.0e0;
    double res = 2.0e0 / phys::math::pi * atan(2 * V[ik] / (V[ii] - V[kk]));
    if (abs(res) < 1.0e-8) res = copysign(1.0e-8, res);
    return res;
}

int calc_wrho(num_complex* wrho,  // distorted rho
              num_complex* rho,   // rho_ele
              double xi,          // xi must be 1
              double gamma,       // gamma must be 0
              double alpha) {
    // initialize distorted-density
    num_real L       = 1.0e0 - log(abs(alpha));  // @NOTE
    num_complex norm = 0.0e0;
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

double calc_Ew(num_real* E, num_complex* wrho, int occ) {
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
int calc_distforce(num_real* f1,       // to be calculated
                   num_real* E,        // (input)
                   num_real* dE,       // (input)
                   num_complex* wrho,  // distorted rho
                   num_complex* rho,   // rho_ele
                   double alpha) {
    // initialize distorted-density
    num_real L            = 1.0e0 - log(abs(alpha));  // @NOTE
    num_real rate_default = (L == 1.0e0) ? 1.0e0 : 0.0e0;

    // averaged adiabatic energy by distorted-density
    double Ew = 0.0e0;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) Ew += std::real(wrho[ii] * E[i]);

    // distortion force (when L= 1 [EHR], the distortion force is zero)
    for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
        num_real* dEj = dE + jFF;
        f1[j]         = 0.0e0;
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

int hopping_impulse(num_real* direction, num_real* np, num_real* nm,  //
                    num_real Efrom, num_real Eto, int from, int to, bool reflect) {
    if (to == from) return from;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    num_real coeffa = 0.0f, coeffb = 0.0f, coeffc = Eto - Efrom;
    for (int i = 0; i < Dimension::N; ++i) {
        coeffa += 0.5f * direction[i] * direction[i] / nm[i];
        coeffb += np[i] / nm[i] * direction[i];
    }
    coeffb /= coeffa, coeffc /= coeffa;  // normalization for safety

    num_real coeffd = coeffb * coeffb - 4 * coeffc;
    if (coeffd > 0) {
        num_real x1 = 0.5f * (-coeffb + sqrt(coeffd)), x2 = 0.5f * (-coeffb - sqrt(coeffd));
        num_real xx = (std::abs(x1) < std::abs(x2)) ? x1 : x2;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return to;
    } else if (reflect) {  // 2008Algorithm
        num_real xx = -coeffb;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return from;
    } else {  // 1990Algorithm, do nothing
        return from;
    }
    return from;
}

void Kernel_Elec_CMSH::read_param_impl(Param* PM) {
    cmsh_type = CMSHPolicy::_from(PM->get<std::string>("cmsh_flag", LOC(), "CVSH"));

    alpha0 = PM->get<num_real>("alpha0", LOC(), 0.5);
    gamma1 = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));

    if (gamma1 < -1.5) gamma1 = Kernel_Elec_MMSH::gamma_opt(Dimension::F);
    if (gamma1 < -0.5) gamma1 = Kernel_Elec_CMM::gamma_wigner(Dimension::F);

    gamma2          = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1             = (1 + Dimension::F * gamma1);
    xi2             = (1 + Dimension::F * gamma2);
    use_wmm         = PM->get<bool>("use_wmm", LOC(), false);
    use_sqc         = PM->get<bool>("use_sqc", LOC(), false);
    use_fssh        = PM->get<bool>("use_fssh", LOC(), false);
    use_strange_win = PM->get<bool>("use_strange_win", LOC(), false);
    use_focus       = PM->get<bool>("use_focus", LOC(), false);
    dt              = PM->get<double>("dt", LOC(), phys::time_d);

    hopping_type1 = 0;
    hopping_type2 = 0;
    reflect       = true;
    use_cv        = false;
    dynamic_alpha = false;
    switch (cmsh_type) {
        case CMSHPolicy::EHR:
            Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;
            break;
        case CMSHPolicy::BOSH:
            Kernel_NADForce::NADForce_type = NADForcePolicy::BO;
            hopping_type1                  = 0;  // max(rho_ele)
            hopping_type2                  = 0;  // MASH direction
            reflect                        = true;
            break;
        case CMSHPolicy::CVSH:
            Kernel_NADForce::NADForce_type = NADForcePolicy::CV;
            hopping_type1                  = 1;  // max(rho_nuc)
            hopping_type2                  = 2;  // P direction
            reflect                        = false;
            use_cv                         = true;
            break;
        case CMSHPolicy::BOSD:
            Kernel_NADForce::NADForce_type = NADForcePolicy::BOSD;
            dynamic_alpha                  = true;
            break;
        case CMSHPolicy::CVSD:
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

void Kernel_Elec_CMSH::init_data_impl(DataSet* DS) {
    alpha     = DS->reg<num_real>("integrator.alpha", Dimension::P);
    Epot      = DS->reg<num_real>("integrator.Epot", Dimension::P);
    p         = DS->reg<num_real>("integrator.p", Dimension::PN);
    m         = DS->reg<num_real>("integrator.m", Dimension::PN);
    fadd      = DS->reg<num_real>("integrator.fadd", Dimension::PN);
    ftmp      = DS->reg<num_real>("integrator.tmp.ftmp", Dimension::N);
    wrho      = DS->reg<num_complex>("integrator.tmp.wrho", Dimension::FF);
    vpes      = DS->reg<num_real>("model.vpes", Dimension::P);
    V         = DS->reg<num_real>("model.V", Dimension::PFF);
    E         = DS->reg<num_real>("model.rep.E", Dimension::PF);
    dE        = DS->reg<num_real>("model.rep.dE", Dimension::PNFF);
    T         = DS->reg<num_real>("model.rep.T", Dimension::PFF);
    H         = DS->reg<num_complex>("model.rep.H", Dimension::PFF);
    direction = DS->reg<num_real>("integrator.tmp.direction", Dimension::N);
}

void Kernel_Elec_CMSH::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* w       = Kernel_Elec::w + iP;
        num_complex* wz_A    = Kernel_Elec::wz_A + iP;
        num_complex* wz_D    = Kernel_Elec::wz_D + iP;
        num_complex* ww_A    = Kernel_Elec::ww_A + iP;
        num_complex* ww_D    = Kernel_Elec::ww_D + iP;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        num_real* T          = Kernel_Elec::T + iP * Dimension::FF;
        int* occ_nuc         = Kernel_Elec::occ_nuc + iP;
        num_real* alpha      = this->alpha + iP;
        num_real* V          = this->V + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////
        alpha[0] = (dynamic_alpha) ? calc_alpha(V) : alpha0;

        w[0]     = num_complex(Dimension::F);  ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;          ///< initial occupation
        if (use_focus) {
            Kernel_Elec_CMM::c_focus(c, xi1, gamma1, Kernel_Elec::occ0, Dimension::F);
        } else if (use_sqc) {
            Kernel_Elec_SQC::c_window(c, Kernel_Elec::occ0, SQCPolicy::TRI,
                                      Dimension::F);  ///< initial c: non-standard c
        } else {
            Kernel_Elec_CMM::c_sphere(c, Dimension::F);  ///< initial c on standard sphere
        }
        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele
        Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, (use_sqc ? 1.0e0 : xi1), gamma1, Dimension::F, use_cv,
                                  *occ_nuc);  ///< initial rho_nuc
        ARRAY_EYE(U, Dimension::F);           ///< initial propagator

        // BO occupation in adiabatic representation
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);
        *occ_nuc = Kernel_Elec_SH::max_choose(rho_nuc);
        if (use_fssh) *occ_nuc = Kernel_Elec_SH::pop_choose(rho_nuc);
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
        int max_pop    = Kernel_Elec_SH::max_choose(rho_ele);
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_A[0]        = 4.0 - 1.0 / (max_val * max_val);
        Kernel_Representation::transform(rho_ele, T, Dimension::F,         //
                                         RepresentationPolicy::Adiabatic,  //
                                         RepresentationPolicy::Diabatic,   //
                                         SpacePolicy::L);
        wz_D[0] = std::abs(rho_ele[0] - rho_ele[3]);
        max_pop = Kernel_Elec_SH::max_choose(rho_ele);
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         RepresentationPolicy::Diabatic,        //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
    }

    _DataSet->set("init.wz_A", Kernel_Elec::wz_A, Dimension::P);
    _DataSet->set("init.wz_D", Kernel_Elec::wz_D, Dimension::P);
    Kernel_Elec::ww_A_init    = _DataSet->set("init.ww_A", Kernel_Elec::ww_A, Dimension::P);
    Kernel_Elec::ww_D_init    = _DataSet->set("init.ww_D", Kernel_Elec::ww_D, Dimension::P);
    Kernel_Elec::c_init       = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->set("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::rho_nuc_init = _DataSet->set("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    Kernel_Elec::T_init       = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);

    exec_kernel(stat);
    // for (int iP = 0; iP < Dimension::P; ++iP) {  // @debug only for scattering problem
    //     num_real* vpes = this->vpes + iP;
    //     num_real* E    = this->E + iP;
    //     num_real* Epot = this->Epot + iP;
    //     num_real* p    = this->p + iP;
    //     num_real* m    = this->m + iP;
    //     double Ekin = 0.0e0;
    //     for (int j = 0; j < Dimension::N; ++j) Ekin += 0.5e0 * p[j] * p[j] / m[j];
    //     double Ekin2 = Ekin + E[Kernel_Elec::occ0] - Epot[0];
    //     double scale = sqrt(std::max({0.0e0, Ekin2 / Ekin}));
    //     for (int j = 0; j < Dimension::N; ++j) p[j] *= scale;
    // }
}

int Kernel_Elec_CMSH::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        int* occ_nuc              = Kernel_Elec::occ_nuc + iP;
        num_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        num_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* rho_nuc_init = Kernel_Elec::rho_nuc_init + iP * Dimension::FF;
        num_real* T               = Kernel_Elec::T + iP * Dimension::FF;
        num_real* T_init          = Kernel_Elec::T_init + iP * Dimension::FF;
        num_complex* K0           = Kernel_Elec::K0 + iP * Dimension::FF;
        num_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        num_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;
        num_complex* K1QA         = Kernel_Elec::K1QA + iP * Dimension::FF;
        num_complex* K2QA         = Kernel_Elec::K2QA + iP * Dimension::FF;
        num_complex* K1DA         = Kernel_Elec::K1DA + iP * Dimension::FF;
        num_complex* K2DA         = Kernel_Elec::K2DA + iP * Dimension::FF;
        num_complex* K1QD         = Kernel_Elec::K1QD + iP * Dimension::FF;
        num_complex* K2QD         = Kernel_Elec::K2QD + iP * Dimension::FF;
        num_complex* K1DD         = Kernel_Elec::K1DD + iP * Dimension::FF;
        num_complex* K2DD         = Kernel_Elec::K2DD + iP * Dimension::FF;
        num_complex* ww_A_init    = Kernel_Elec::ww_A_init + iP;
        num_complex* ww_D_init    = Kernel_Elec::ww_D_init + iP;
        num_complex* ww_A         = Kernel_Elec::ww_A + iP;
        num_complex* ww_D         = Kernel_Elec::ww_D + iP;

        num_real* alpha = this->alpha + iP;
        num_real* Epot  = this->Epot + iP;
        num_real* vpes  = this->vpes + iP;
        num_real* V     = this->V + iP * Dimension::FF;
        num_real* E     = this->E + iP * Dimension::F;
        num_real* dE    = this->dE + iP * Dimension::NFF;
        num_real* p     = this->p + iP * Dimension::N;
        num_real* m     = this->m + iP * Dimension::N;
        num_real* fadd  = this->fadd + iP * Dimension::N;
        num_complex* H  = this->H + iP * Dimension::FF;

        //////////////////////////////////////////////////////////////////////

        // * additional evolution can be appended here

        // 1) transform from inp_repr => ele_repr
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = rho_nuc_init[ik];
        Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(rho_nuc, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);

        // 2) propagte along ele_repr
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

        // step 1: determine where to hop (BOSH & CVSH)
        /// 1.1 calc Efrom
        num_real Efrom, Eto;
        Efrom = calc_Ew(E, rho_nuc, *occ_nuc);
        /// 1.2 calc to
        int to;
        switch (hopping_type1) {
            case 0: {
                to = Kernel_Elec_SH::max_choose(rho_ele);
                break;
            }
            case 1: {
                to = Kernel_Elec_SH::max_choose(rho_nuc);
                break;
            }
            case 2: {
                to = Kernel_Elec_SH::pop_choose(rho_ele);
                break;
            }
            case 3: {
                to = Kernel_Elec_SH::hopping_choose(rho_ele, H, *occ_nuc, dt);
                break;
            }
            case 4: {
                to = Kernel_Elec_SH::pop_choose(rho_nuc);
                break;
            }
            case 5: {
                to = Kernel_Elec_SH::pop_neg_choose(rho_nuc);
                break;
            }
        }
        Eto = calc_Ew(E, rho_nuc, to);
        // step 2: determine direction to hop
        switch (hopping_type2) {
            case 0: {
                Kernel_Elec_MMSH::hopping_direction(direction, E, dE, rho_ele, *occ_nuc, to);
                break;
            }
            case 1: {
                Kernel_Elec_SH::hopping_direction(direction, dE, *occ_nuc, to);
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

        // smooth dynamics (BOSD & CVSD)
        if (cmsh_type == CMSHPolicy::BOSD || cmsh_type == CMSHPolicy::CVSD) {
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
                double xsolve = (Ew_new - Ew_old) / dt / vdotd;
                for (int j = 0; j < Dimension::N; ++j) fadd[j] += xsolve * dE[j * Dimension::FF + 1];
            }
            for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = wrho[ik];
        }

        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        // 4) calculated TCF in adiabatic rep & diabatic rep respectively
        // 4-1) Adiabatic rep
        int max_pop    = Kernel_Elec_SH::max_choose(rho_ele);
        double max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_A[0]        = 4.0 - 1.0 / (max_val * max_val);
        ww_A[0]        = std::min({abs(ww_A[0]), abs(ww_A_init[0])});

        Kernel_Elec::ker_from_rho(K1QA, rho_ele, 1, 0, Dimension::F, true, max_pop);
        Kernel_Elec::ker_from_rho(K2QA, rho_ele, 1, 0, Dimension::F, true, max_pop);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QA[i * Dimension::Fadd1] = (abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        if (use_fssh) { Kernel_Elec::ker_from_rho(K2QA, rho_ele, 1, 0, Dimension::F, true, *occ_nuc); }

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

        max_pop = Kernel_Elec_SH::max_choose(rho_ele);  // (in dia rep)
        max_val = std::abs(rho_ele[max_pop * Dimension::Fadd1]);
        ww_D[0] = 4.0 - 1.0 / (max_val * max_val);
        ww_D[0] = std::min({abs(ww_D[0]), abs(ww_D_init[0])});

        Kernel_Elec::ker_from_rho(K1QD, rho_ele, 1, 0, Dimension::F, true, max_pop);
        Kernel_Elec::ker_from_rho(K2QD, rho_ele, 1, 0, Dimension::F);
        for (int i = 0; i < Dimension::F; ++i) {
            K2QD[i * Dimension::Fadd1] = (abs(rho_ele[i * Dimension::Fadd1]) < 1 / xi1) ? 0.0e0 : 1.0e0;
        }
        if (use_strange_win) calc_wrho(K2QD, rho_ele, 1, 0, 0.2);
        if (use_sqc) { Kernel_Elec_SQC::ker_binning(K2QD, rho_ele, SQCPolicy::TRI); }

        ARRAY_MAT_DIAG(K1DD, K1QD, Dimension::F);
        ARRAY_MAT_DIAG(K2DD, K2QD, Dimension::F);

        // 5) transform back from tcf_repr => inp_repr
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

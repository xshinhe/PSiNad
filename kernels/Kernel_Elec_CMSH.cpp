#include "Kernel_Elec_CMSH.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec_CMM.h"
#include "Kernel_Elec_MMSH.h"
#include "Kernel_Elec_SH.h"
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

int Kernel_Elec_CMSH::hopping_impulse(num_real* direction, num_real* np, num_real* nm,  //
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
    gamma1  = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));
    gamma2  = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1     = (1 + Dimension::F * gamma1);
    xi2     = (1 + Dimension::F * gamma2);
    use_cv  = PM->get<bool>("use_cv", LOC(), true);
    use_wmm = PM->get<bool>("use_wmm", LOC(), false);
    reflect = PM->get<bool>("reflect", LOC(), true);  ///< reflect scheme in hopping dynamics

    dt            = PM->get<double>("dt", LOC(), phys::time_d);
    hopping_type1 = PM->get<int>("hopping_type1", LOC(), 0);
    hopping_type2 = PM->get<int>("hopping_type2", LOC(), 0);
    hopping_type3 = PM->get<int>("hopping_type3", LOC(), 0);
}

void Kernel_Elec_CMSH::init_data_impl(DataSet* DS) {
    p         = DS->reg<num_real>("integrator.p", Dimension::PN);
    m         = DS->reg<num_real>("integrator.m", Dimension::PN);
    E         = DS->reg<num_real>("model.rep.E", Dimension::PF);
    dE        = DS->reg<num_real>("model.rep.dE", Dimension::PNFF);
    T         = DS->reg<num_real>("model.rep.T", Dimension::PFF);
    H         = DS->reg<num_complex>("model.rep.H", Dimension::PFF);
    direction = DS->reg<num_real>("integrator.tmp.direction", Dimension::N);
}

void Kernel_Elec_CMSH::init_calc_impl(int stat) {
    Kernel_NADForce::NADForce_type = NADForcePolicy::BO;

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

        /////////////////////////////////////////////////////////////////

        w[0]     = num_complex(Dimension::F);                     ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;                             ///< initial occupation
        Kernel_Elec_CMM::c_sphere(c, Dimension::F);               ///< initial c on standard sphere
        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele
        Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, *occ_nuc);  ///< initial rho_nuc
        ARRAY_EYE(U, Dimension::F);  ///< initial propagator

        // BO occupation in adiabatic representation
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::inp_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);
        *occ_nuc = Kernel_Elec_SH::max_choose(rho_nuc);
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

        num_real* E    = this->E + iP * Dimension::F;
        num_real* dE   = this->dE + iP * Dimension::NFF;
        num_real* p    = this->p + iP * Dimension::N;
        num_real* m    = this->m + iP * Dimension::N;
        num_complex* H = this->H + iP * Dimension::FF;

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

        // step 1: determine where to hop
        num_real Efrom, Eto;
        switch (Kernel_NADForce::NADForce_type) {
            case NADForcePolicy::BO: {
                Efrom = E[*occ_nuc];
                break;
            }
            default: {
                Efrom = std::real(ARRAY_TRACE2(rho_nuc, H, Dimension::F, Dimension::F));
                break;
            }
        }
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
        }
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
        switch (hopping_type3) {
            case 0: {  // always BO
                Kernel_NADForce::NADForce_type = NADForcePolicy::BO;
                Eto                            = E[to];
                break;
            }
            case 1: {  // always EHR
                Kernel_NADForce::NADForce_type = NADForcePolicy::EHR;
                Kernel_Elec::ker_from_rho(rho_nuc, rho_ele_init, xi1, gamma1, Dimension::F, use_cv,
                                          to);                                          ///< re-initial rho_nuc
                Kernel_Representation::transform(rho_nuc, T_init, Dimension::F,         //
                                                 Kernel_Representation::inp_repr_type,  //
                                                 Kernel_Representation::ele_repr_type,  //
                                                 SpacePolicy::L);
                ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                                 Kernel_Representation::ele_repr_type,  //
                                                 RepresentationPolicy::Adiabatic,       //
                                                 SpacePolicy::L);
                Eto = std::real(ARRAY_TRACE2(rho_nuc, H, Dimension::F, Dimension::F));
                break;
            }
        }
        *occ_nuc = hopping_impulse(direction, p, m, Efrom, Eto, *occ_nuc, to, reflect);
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
        if (abs(rho_ele[max_pop * Dimension::Fadd1]) < 1) K2QA[max_pop * Dimension::Fadd1] = 0.0e0;

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
        Kernel_Elec::ker_from_rho(K2QD, rho_ele, 1, 0, Dimension::F, true, max_pop);
        if (abs(rho_ele[max_pop * Dimension::Fadd1]) < 1 / xi1) K2QD[max_pop * Dimension::Fadd1] = 0.0e0;

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

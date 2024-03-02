#include "Kernel_Elec_SH.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
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

namespace kids {

int Kernel_Elec_SH::max_choose(kids_complex* rho) {
    int imax       = 0;
    kids_real vmax = 0.0f;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        if (std::real(rho[ii]) > vmax) {
            vmax = std::real(rho[ii]);
            imax = i;
        }
    }
    return imax;
}

int Kernel_Elec_SH::pop_choose(kids_complex* rho) {
    kids_real rand_tmp;
    kids_real sum = 0.0f;
    Kernel_Random::rand_uniform(&rand_tmp);
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        sum += std::min({std::max({std::real(rho[ii]), 0.0e0}), 1.0e0});
        if (rand_tmp < sum) return i;
    }
    return 0;
}

int Kernel_Elec_SH::pop_neg_choose(kids_complex* rho) {
    kids_real rand_tmp;
    kids_real sum = 0.0f;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) sum += std::abs(rho[ii]);
    Kernel_Random::rand_uniform(&rand_tmp);
    rand_tmp *= sum;
    sum = 0.0e0;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        sum += std::abs(rho[ii]);
        if (rand_tmp < sum) return i;
    }
    return 0;
}

int Kernel_Elec_SH::hopping_choose(kids_complex* rho, kids_complex* H, int from, kids_real dt) {
    int to = from;
    kids_real rand_tmp, sumprob = 0.0f;
    kids_real rhoii = std::real(rho[from * Dimension::Fadd1]);
    Kernel_Random::rand_uniform(&rand_tmp);

    for (int n = 0; n < Dimension::F; ++n) {
        kids_real prob =
            (n == from) ? 0.0f
                        : -2.0f * std::imag(rho[n * Dimension::F + from] * H[from * Dimension::F + n]) / rhoii * dt;
        prob = (prob > 1.0f) ? 1.0f : ((prob < 0.0f) ? 0.0f : prob);  // hopping cut-off
        sumprob += prob;
        if (rand_tmp < sumprob) {
            to = n;
            break;
        }
    }
    return to;
}

void Kernel_Elec_SH::hopping_direction(kids_real* direction, kids_real* dE, int from, int to) {
    if (to == from) return;
    for (int i = 0; i < Dimension::N; ++i) { direction[i] = dE[i * Dimension::FF + from * Dimension::F + to]; }
}

int Kernel_Elec_SH::hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm, kids_real* E,  //
                                    int from, int to, bool reflect) {
    if (to == from) return from;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    kids_real coeffa = 0.0f, coeffb = 0.0f, coeffc = E[to] - E[from];
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


void Kernel_Elec_SH::read_param_impl(Param* PM) {
    sh_type = SHPolicy::_from(PM->get<std::string>("sh_flag", LOC(), "FSSH"));
    reflect = PM->get<bool>("reflect", LOC(), false);
    dt      = PM->get<double>("dt", LOC(), phys::time_d);

    std::string REP_CHECK = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    assert(REP_CHECK == "Adiabatic");
}

void Kernel_Elec_SH::init_data_impl(DataSet* DS) {
    p  = DS->def<kids_real>("integrator.p", Dimension::PN);
    m  = DS->def<kids_real>("integrator.m", Dimension::PN);
    E  = DS->def<kids_real>("model.rep.E", Dimension::PF);
    dE = DS->def<kids_real>("model.rep.dE", Dimension::PNFF);
    T  = DS->def<kids_real>("model.rep.T", Dimension::PFF);
    H  = DS->def<kids_complex>("model.rep.H", Dimension::PFF);

    direction = DS->def<kids_real>("integrator.tmp.direction", Dimension::N);
}

void Kernel_Elec_SH::init_calc_impl(int stat) {
    Kernel_NADForce::NADForce_type = NADForcePolicy::BO;  // set BO dynamics for SH

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_complex* w       = Kernel_Elec::w + iP;
        kids_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        kids_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        kids_complex* K1      = Kernel_Elec::K1 + iP * Dimension::FF;
        kids_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc          = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0]     = 1.0e0;                                                               ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;                                                   ///< initial occupation
        for (int i = 0; i < Dimension::F; ++i) c[i] = (i == *occ_nuc) ? 1.0e0 : 0.0e0;  ///< initial c
        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);                        ///< initial rho_ele
        rho_nuc;                                                                        /// not used
        ARRAY_EYE(U, Dimension::F);                                                     ///< initial propagator

        if (Kernel_Representation::inp_repr_type == RepresentationPolicy::Diabatic) {
            ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
            *occ_nuc = pop_choose(Kernel_Elec::rho_ele);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
        }
    }
    Kernel_Elec::c_init       = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->set("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::T_init       = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);
    exec_kernel(stat);
}

int Kernel_Elec_SH::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        int* occ_nuc               = Kernel_Elec::occ_nuc + iP;
        kids_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        kids_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        kids_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        kids_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;
        kids_real* T               = Kernel_Elec::T + iP * Dimension::FF;
        kids_real* T_init          = Kernel_Elec::T_init + iP * Dimension::FF;


        kids_real* E    = this->E + iP * Dimension::F;
        kids_real* dE   = this->dE + iP * Dimension::NFF;
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* m    = this->m + iP * Dimension::N;
        kids_complex* H = this->H + iP * Dimension::FF;

        //////////////////////////////////////////////////////////////////////

        // * additional evolution can be appended here

        for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];
        // 1) transform from inp_repr => ele_repr
        Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);

        // 2) propagte along ele_repr
        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);

        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         RepresentationPolicy::Adiabatic,       //
                                         SpacePolicy::L);

        int to = hopping_choose(rho_ele, H, *occ_nuc, dt);                      // step 1: determine where to hop
        hopping_direction(direction, dE, *occ_nuc, to);                         // step 2: determine direction to hop
        *occ_nuc = hopping_impulse(direction, p, m, E, *occ_nuc, to, reflect);  // step 3: try hop
        Kernel_Elec::ker_from_rho(K1, rho_ele, 1, 0, Dimension::F);
        Kernel_Elec::ker_from_rho(K2, rho_ele, 1, 0, Dimension::F, true, *occ_nuc);


        // 3) transform back from ele_repr => inp_repr
        Kernel_Representation::transform(K1, T, Dimension::F,                   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(K2, T, Dimension::F,                   //
                                         RepresentationPolicy::Adiabatic,       //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
    }
    return stat;
}


};  // namespace kids

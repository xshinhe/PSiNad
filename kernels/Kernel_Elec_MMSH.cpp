#include "Kernel_Elec_MMSH.h"

#include "Kernel_Declare.h"
#include "Kernel_Elec_CMM.h"
#include "Kernel_Elec_SH.h"
#include "Kernel_NADForce.h"
#include "Kernel_Random.h"

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


void Kernel_Elec_MMSH::hopping_direction(num_real* direction, num_real* E, num_real* dE, num_complex* rho, int from,
                                         int to) {
    if (to == from) return;

    for (int i = 0; i < Dimension::N; ++i) {
        direction[i] = 0.0f;
        for (int k = 0; k < Dimension::F; ++k)
            if (k != from)
                direction[i] +=
                    std::real(rho[from * Dimension::F + k] * dE[i * Dimension::FF + k * Dimension::F + from]) /
                    (E[from] - E[k]);
        for (int k = 0; k < Dimension::F; ++k)
            if (k != to)
                direction[i] -= std::real(rho[to * Dimension::F + k] * dE[i * Dimension::FF + k * Dimension::F + to]) /
                                (E[to] - E[k]);
    }
}

void Kernel_Elec_MMSH::read_param_impl(Param* PM) {
    mmsh_type = MMSHPolicy::_from(PM->get<std::string>("mmsh_flag", LOC(), "MASH1"));
    sumover   = PM->get<bool>("sumover", LOC(), true);   ///< sum over sub-manifold Mi other than single Mocc
    focused   = PM->get<bool>("focused", LOC(), false);  ///< focus on sub-manifold Mi
    hopping   = PM->get<bool>("hopping", LOC(), true);   ///< do hopping dynamics
    reflect   = PM->get<bool>("reflect", LOC(), true);   ///< reflect scheme in hopping dynamics

    dt = PM->get<double>("dt", LOC(), phys::time_d);

    std::string REP_CHECK = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    assert(REP_CHECK == "Adiabatic");

    gamma = PM->get<double>("gamma0", LOC(), gamma_opt(Dimension::F));
    xi    = (1 + Dimension::F * gamma);
}

void Kernel_Elec_MMSH::init_data_impl(DataSet* DS) {
    x         = DS->reg<num_real>("integrator.x", Dimension::N);
    p         = DS->reg<num_real>("integrator.p", Dimension::N);
    m         = DS->reg<num_real>("integrator.m", Dimension::N);
    direction = DS->reg<num_real>("integrator.direction", Dimension::N);
    E         = DS->reg<num_real>("model.rep.E", Dimension::F);
    dE        = DS->reg<num_real>("model.rep.dE", Dimension::NFF);
    T         = DS->reg<num_real>("model.rep.T", Dimension::FF);
    H         = DS->reg<num_complex>("model.rep.H", Dimension::FF);

    w_CC = DS->reg<num_complex>("integrator.w_CC");
    w_CP = DS->reg<num_complex>("integrator.w_CP");
    w_PP = DS->reg<num_complex>("integrator.w_PP");
    w_AA = DS->reg<num_complex>("integrator.w_AA");
    w_AD = DS->reg<num_complex>("integrator.w_AD");
    w_DD = DS->reg<num_complex>("integrator.w_DD");

    DS->set("init.w_CC", w_CC);
    DS->set("init.w_CP", w_CP);
    DS->set("init.w_PP", w_PP);
    DS->set("init.w_AA", w_AA);
    DS->set("init.w_AD", w_AD);
    DS->set("init.w_DD", w_DD);
}

void Kernel_Elec_MMSH::init_calc_impl(int stat) {
    /**
     * 1) if do hopping dynamics, here adopt the strategy of MASH
     * 2) otherwise, do mean-force dynamics acccording to
     */
    Kernel_NADForce::NADForce_type = (hopping) ? NADForcePolicy::BO : NADForcePolicy::EHR;

    // sampling electonic DOFs in <ini_repr>
    switch (focused) {
        case true: {
            double randu;
            do {
                Kernel_Random::rand_uniform(&randu);
                *Kernel_Elec::occ_nuc = int(randu * Dimension::F);
            } while (!sumover && (*Kernel_Elec::occ_nuc) != (Kernel_Elec::occ0));

            for (int i = 0; i < Dimension::F; ++i) {
                Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
                Kernel_Elec::c[i] = (i == *Kernel_Elec::occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) *
                                    (cos(randu) + phys::math::im * sin(randu));
            }
            Kernel_Elec::ker_from_c(Kernel_Elec::rho_ele, Kernel_Elec::c, xi, gamma, Dimension::F);
            break;
        }
        case false: {
            do {
                Kernel_Elec_CMM::c_sphere(Kernel_Elec::c, Dimension::F);
                Kernel_Elec::ker_from_c(Kernel_Elec::rho_ele, Kernel_Elec::c, 1, 0, Dimension::F);
                *Kernel_Elec::occ_nuc = Kernel_Elec_SH::max_choose(Kernel_Elec::rho_ele);
            } while (!sumover && (*Kernel_Elec::occ_nuc) != (Kernel_Elec::occ0));
            break;
        }
    }
    // ARRAY_SHOW(Kernel_Elec::rho_ele, Dimension::F, Dimension::F);

    // REP: <ini_repr> ==> ADIABATIC
    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
    }

    // ARRAY_SHOW(Kernel_Elec::rho_ele, Dimension::F, Dimension::F);

    // determine the weighting factor
    Kernel_Elec::w[0] = (sumover) ? num_complex(Dimension::F) : 1.0e0;
    switch (mmsh_type) {
        case MMSHPolicy::MASH1:
            w_CC[0] = Kernel_Elec::w[0] * 3.0e0;
            w_CP[0] = Kernel_Elec::w[0] * 2.0e0;
            w_PP[0] = Kernel_Elec::w[0] * 2.0e0 * std::abs(Kernel_Elec::rho_ele[0] - Kernel_Elec::rho_ele[3]);
            break;
        case MMSHPolicy::MASH2:
            // w_CC = Kernel_Elec::w[0] * 3.0e0;
            // w_CP = Kernel_Elec::w[0] * (1 + Dimension::F * gamma);
            // w_PP = Kernel_Elec::w[0];
            w_CC[0] = Kernel_Elec::w[0];
            w_CP[0] = Kernel_Elec::w[0];
            w_PP[0] = Kernel_Elec::w[0];
            break;
    }
    w_AA[0] = w_CC[0];
    w_AD[0] = w_CP[0] - w_CC[0];
    w_DD[0] = w_PP[0] - 2.0e0 * w_CP[0] + w_CC[0];
    _DataSet->set("init.w_CC", w_CC, 1);
    _DataSet->set("init.w_CP", w_CP, 1);
    _DataSet->set("init.w_PP", w_PP, 1);
    _DataSet->set("init.w_AA", w_AA, 1);
    _DataSet->set("init.w_AD", w_AD, 1);
    _DataSet->set("init.w_DD", w_DD, 1);

    *Kernel_Elec::occ_nuc = Kernel_Elec_SH::max_choose(Kernel_Elec::rho_ele);

    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::K1, T, Kernel_Elec::K1, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::K1Q, T, Kernel_Elec::K1Q, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
    }

    for (int ik = 0; ik < Dimension::FF; ++ik) Kernel_Elec::rho_nuc[ik] = Kernel_Elec::rho_ele[ik];
    exec_kernel(stat);
}

int Kernel_Elec_MMSH::exec_kernel_impl(int stat) {
    /**
     * additional evolution appended to Kernel_Update_rho (auxilrary variables)
     */

    // REP: <ini_repr> ==> ADIABATIC
    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
    }

    // step 1: determine where to hop
    int to = Kernel_Elec_SH::max_choose(Kernel_Elec::rho_ele);
    // step 2: determine direction to hop
    hopping_direction(direction, E, dE, Kernel_Elec::rho_ele, *Kernel_Elec::occ_nuc, to);
    // step 3: try hop
    *Kernel_Elec::occ_nuc = Kernel_Elec_SH::hopping_impulse(direction, p, m, E, *Kernel_Elec::occ_nuc, to, reflect);

    // K1 & K1Q
    Kernel_Elec::ker_from_rho_quantize(Kernel_Elec::K1, Kernel_Elec::rho_ele, 1, 0, *Kernel_Elec::occ_nuc,
                                       Dimension::F);
    ARRAY_CLEAR(Kernel_Elec::K1Q, Dimension::FF);
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) Kernel_Elec::K1Q[ii] = Kernel_Elec::K1[ii];

    // K2
    Kernel_Elec::ker_from_rho(Kernel_Elec::K2, Kernel_Elec::rho_ele, xi, gamma, Dimension::F);

    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::K1, T, Kernel_Elec::K1, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::K1Q, T, Kernel_Elec::K1Q, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::K2, T, Kernel_Elec::K2, T, Dimension::F, Dimension::F, Dimension::F,
                             Dimension::F);
    }

    int imax_init_rep = Kernel_Elec_SH::max_choose(Kernel_Elec::rho_ele);
    Kernel_Elec::ker_from_rho_quantize(Kernel_Elec::K2Q, Kernel_Elec::rho_ele, xi, gamma, imax_init_rep, Dimension::F);

    if (!hopping) {
        for (int ik = 0; ik < Dimension::FF; ++ik) {
            Kernel_Elec::rho_nuc[ik] = Kernel_Elec::K2[ik];  //
            // Kernel_Elec::rho_nuc[ik] = Kernel_Elec::K2Q[ik];  //
        }
    }

    return stat;
}


};  // namespace PROJECT_NS

#include "Kernel_Elec_MMSH.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec_CMM.h"
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


void Kernel_Elec_MMSH::hopping_direction(kids_real* direction, kids_real* E, kids_real* dE, kids_complex* rho, int from,
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
    use_cv    = PM->get<bool>("use_cv", LOC(), false);   ///< use cv if not hopping

    dt = PM->get<double>("dt", LOC(), phys::time_d);

    std::string REP_CHECK = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    assert(REP_CHECK == "Adiabatic");

    gamma = PM->get<double>("gamma0", LOC(), gamma_opt(Dimension::F));
    xi    = (1 + Dimension::F * gamma);
}

void Kernel_Elec_MMSH::init_data_impl(DataSet* DS) {
    p  = DS->def<kids_real>("integrator.p", Dimension::PN);
    m  = DS->def<kids_real>("integrator.m", Dimension::PN);
    E  = DS->def<kids_real>("model.rep.E", Dimension::PF);
    dE = DS->def<kids_real>("model.rep.dE", Dimension::PNFF);
    T  = DS->def<kids_real>("model.rep.T", Dimension::PFF);
    H  = DS->def<kids_complex>("model.rep.H", Dimension::PFF);

    direction = DS->def<kids_real>("integrator.tmp.direction", Dimension::N);
}

void Kernel_Elec_MMSH::init_calc_impl(int stat) {
    Kernel_NADForce::NADForce_type = (hopping) ? NADForcePolicy::BO : NADForcePolicy::EHR;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_complex* w       = Kernel_Elec::w + iP;
        kids_complex* w_CC    = Kernel_Elec::w_CC + iP;
        kids_complex* w_CP    = Kernel_Elec::w_CP + iP;
        kids_complex* w_PP    = Kernel_Elec::w_PP + iP;
        kids_complex* w_AA    = Kernel_Elec::w_AA + iP;
        kids_complex* w_AD    = Kernel_Elec::w_AD + iP;
        kids_complex* w_DD    = Kernel_Elec::w_DD + iP;
        kids_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        kids_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        kids_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc          = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0] = (sumover) ? kids_complex(Dimension::F) : 1.0e0;  ///< initial measure

        do {  // sampling occ_nuc, c, rho_ele [branch for sumover & focused]
            if (focused) {
                double randu;
                Kernel_Random::rand_uniform(&randu);
                *occ_nuc = int(randu * Dimension::F);
                for (int i = 0; i < Dimension::F; ++i) {
                    Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
                    c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu);
                }
                Kernel_Elec::ker_from_c(rho_ele, c, xi, gamma, Dimension::F);
            } else {
                Kernel_Elec_CMM::c_sphere(c, Dimension::F);
                Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);
                *occ_nuc = Kernel_Elec_SH::max_choose(rho_ele);
            }
        } while (!sumover && (*occ_nuc) != (Kernel_Elec::occ0));

        Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi, gamma, Dimension::F, use_cv, *occ_nuc);  ///< initial rho_nuc
        ARRAY_EYE(U, Dimension::F);                                                              ///< initial propagator

        // update initial weighting factor & choose occupied state in adiabatic rep
        if (Kernel_Representation::inp_repr_type == RepresentationPolicy::Diabatic) {
            ARRAY_MATMUL3_TRANS1(rho_ele, T, rho_ele, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        }
        switch (mmsh_type) {
            case MMSHPolicy::MASH1:
                w_CC[0] = w[0] * 3.0e0;
                w_CP[0] = w[0] * 2.0e0;
                w_PP[0] = w[0] * 2.0e0 * std::abs(rho_ele[0] - rho_ele[3]);
                break;
            case MMSHPolicy::MASH2:
                // w_CC = w[0] * 3.0e0;
                // w_CP = w[0] * (1 + Dimension::F * gamma);
                // w_PP = w[0];
                w_CC[0] = w[0];
                w_CP[0] = w[0];
                w_PP[0] = w[0];
                break;
        }
        w_AA[0] = w_CC[0];
        w_AD[0] = w_CP[0] - w_CC[0];
        w_DD[0] = w_PP[0] - 2.0e0 * w_CP[0] + w_CC[0];

        *occ_nuc = Kernel_Elec_SH::max_choose(rho_ele);

        if (Kernel_Representation::inp_repr_type == RepresentationPolicy::Diabatic) {
            ARRAY_MATMUL3_TRANS2(rho_ele, T, rho_ele, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        }
    }

    _DataSet->def("init.w_CC", Kernel_Elec::w_CC, Dimension::P);
    _DataSet->def("init.w_CP", Kernel_Elec::w_CP, Dimension::P);
    _DataSet->def("init.w_PP", Kernel_Elec::w_PP, Dimension::P);
    _DataSet->def("init.w_AA", Kernel_Elec::w_AA, Dimension::P);
    _DataSet->def("init.w_AD", Kernel_Elec::w_AD, Dimension::P);
    _DataSet->def("init.w_DD", Kernel_Elec::w_DD, Dimension::P);
    Kernel_Elec::c_init       = _DataSet->def("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->def("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::rho_nuc_init = _DataSet->def("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    exec_kernel(stat);
}

int Kernel_Elec_MMSH::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        int* occ_nuc               = Kernel_Elec::occ_nuc + iP;
        kids_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        kids_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        kids_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        kids_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        kids_complex* rho_nuc_init = Kernel_Elec::rho_nuc_init + iP * Dimension::FF;
        kids_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        kids_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;
        kids_complex* K1DA         = Kernel_Elec::K1DA + iP * Dimension::FF;
        kids_complex* K2DA         = Kernel_Elec::K2DA + iP * Dimension::FF;

        kids_real* E    = this->E + iP * Dimension::F;
        kids_real* T    = this->T + iP * Dimension::FF;
        kids_real* dE   = this->dE + iP * Dimension::NFF;
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* m    = this->m + iP * Dimension::N;
        kids_complex* H = this->H + iP * Dimension::FF;

        //////////////////////////////////////////////////////////////////////

        // * additional evolution can be appended here

        if (Kernel_Representation::inp_repr_type == RepresentationPolicy::Diabatic) {
            ARRAY_MATMUL_TRANS1(U, T, U, Dimension::F, Dimension::F, Dimension::F);
        }

        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele_init, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc_init, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);

        int to = Kernel_Elec_SH::max_choose(rho_ele);                // step 1: determine where to hop
        hopping_direction(direction, E, dE, rho_ele, *occ_nuc, to);  // step 2: determine direction to hop
        *occ_nuc = Kernel_Elec_SH::hopping_impulse(direction, p, m, E, *occ_nuc, to, reflect);  // step 3: try hop

        Kernel_Elec::ker_from_rho(K1, rho_ele, 1, 0, Dimension::F, true, *occ_nuc);
        ARRAY_MAT_DIAG(K1DA, K1, Dimension::F);
        Kernel_Elec::ker_from_rho(K2, rho_ele, xi, gamma, Dimension::F);

        if (Kernel_Representation::inp_repr_type == RepresentationPolicy::Diabatic) {
            ARRAY_MATMUL(U, T, U, Dimension::F, Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(rho_ele, T, rho_ele, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(K1, T, K1, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(K1DA, T, K1DA, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(K2, T, K2, T, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        }
        int imax_init_rep = Kernel_Elec_SH::max_choose(rho_ele);
        Kernel_Elec::ker_from_rho(K2DA, rho_ele, xi, gamma, Dimension::F, true, imax_init_rep);
    }
    return stat;
}


};  // namespace PROJECT_NS

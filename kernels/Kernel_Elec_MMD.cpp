#include "Kernel_Elec_MMD.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
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

/**
 * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
 */
int Kernel_Elec_MMD::rho_focus(num_complex* rho, int iocc, double gamma_ou, double gamma_uu, int fdim, bool rand_act,
                               bool pure_phase, bool cont_phase) {
    double randu = 1.0e0;
    if (pure_phase) {
        for (int j = 0; j < fdim; ++j) {
            if (j == iocc) {
                rho[iocc * fdim + j] = 1.0e0;
                continue;
            }
            Kernel_Random::rand_uniform(&randu);
            randu = (cont_phase) ? phys::math::twopi * randu : phys::math::halfpi * (int(randu / 0.25f) + 1);
            rho[iocc * fdim + j] = cos(randu) + phys::math::im * sin(randu);
            rho[j * fdim + iocc] = std::conj(rho[iocc * fdim + j]);
        }
        for (int i = 0, ij = 0; i < fdim; ++i) {
            for (int j = 0; j < fdim; ++j, ++ij) {
                if (i == iocc || j == iocc) continue;
                rho[ij] = rho[iocc * fdim + j] / rho[iocc * fdim + i];
            }
        }
    } else {
        for (int i = 0; i < fdim; ++i) {
            for (int j = i + 1; j < fdim; ++j) {
                Kernel_Random::rand_uniform(&randu);
                randu = (cont_phase) ? phys::math::twopi * randu : phys::math::halfpi * (int(randu / 0.25f) + 1);
                rho[i * fdim + j] = cos(randu) + phys::math::im * sin(randu);
                rho[j * fdim + i] = std::conj(rho[i * fdim + j]);
            }
        }
    }
    Kernel_Random::rand_uniform(&randu);
    double occrand = (rand_act) ? int(randu * fdim) : iocc;
    for (int i = 0, ij = 0; i < fdim; ++i) {
        for (int j = 0; j < fdim; ++j, ++ij) {
            if (i == j) {
                rho[ij] = (i == occrand) ? phys::math::iu : phys::math::iz;
            } else if (i == occrand || j == occrand) {
                rho[ij] *= gamma_ou;
            } else {
                rho[ij] *= gamma_uu;
            }
        }
    }
    return 0;
}

void Kernel_Elec_MMD::read_param_impl(Param* PM) {
    mmd_type = MMDPolicy::_from(PM->get<std::string>("mmd_flag", LOC(), "MMF"));
    scale    = PM->get<num_real>("scale", LOC(), 1.0);
    Fref     = PM->get<int>("Fref", LOC(), Dimension::F);
    switch (mmd_type) {
        case MMDPolicy::MMF:
            gamma_uu = (std::sqrt(scale * Fref + 1) - 1) / Fref;
            gamma_ou = sqrt(gamma_uu * (1.0f + gamma_uu));
            break;
        case MMDPolicy::TWA:
            gamma_uu = 0.0f;
            gamma_ou = sqrt(0.5f * scale);
            break;
        case MMDPolicy::MID: {
            double scale_k = PM->get<num_real>("scale_k", LOC(), 1.0);
            double r       = (Fref + 2 * scale_k) / (1 + scale_k);
            gamma_uu       = (std::sqrt(r / (1 + scale_k) * scale + 1) - 1) / r;
            gamma_ou       = std::sqrt((1 + scale_k) * gamma_uu * (1 + gamma_uu));
            break;
        }
    }
    double gamma_uu_in = PM->get<num_real>("gamma", LOC(), -1.0);
    if (gamma_uu_in >= 0) {
        gamma_uu           = gamma_uu_in;
        double gamma_ou_in = PM->get<num_real>("gamma_ou", LOC(), -1.0);
        gamma_ou           = (gamma_uu_in >= 0) ? gamma_uu_in : sqrt((1 + gamma_uu) * gamma_uu);
    }

    pure_phase = PM->get<bool>("pure_phase", LOC(), true);
    cont_phase = PM->get<bool>("cont_phase", LOC(), true);
    rand_act   = PM->get<bool>("rand_act", LOC(), false);
}

void Kernel_Elec_MMD::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* w       = Kernel_Elec::w + iP;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc         = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0]     = (rand_act) ? num_complex(Dimension::F) : 1.0e0;  ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;                               ///< initial occupation
        c;
        rho_focus(rho_ele, *occ_nuc, gamma_ou, gamma_uu, Dimension::F, rand_act, pure_phase, cont_phase);
        rho_nuc;                     ///< initial rho_nuc (not used)
        ARRAY_EYE(U, Dimension::F);  ///< initial propagator
    }
    Kernel_Elec::c_init       = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->set("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    // Kernel_Elec::rho_nuc_init = _DataSet->set("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    Kernel_Elec::T_init = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);
    exec_kernel(stat);
}

int Kernel_Elec_MMD::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* c            = Kernel_Elec::c + iP * Dimension::F;
        num_complex* c_init       = Kernel_Elec::c_init + iP * Dimension::F;
        num_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        num_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        num_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;

        num_real* T      = Kernel_Elec::T + iP * Dimension::FF;
        num_real* T_init = Kernel_Elec::T_init + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////

        for (int i = 0; i < Dimension::F; ++i) c[i] = c_init[i];
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];

        // 1) transform from inp_repr => ele_repr
        Kernel_Representation::transform(c, T_init, Dimension::F,               //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::H);
        Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);
        // 2) propagte along ele_repr
        ARRAY_MATMUL(c, U, c, Dimension::F, Dimension::F, 1);
        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);

        // 3) transform back from ele_repr => inp_repr
        Kernel_Representation::transform(c, T, Dimension::F,                    //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::H);
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = rho_ele[ik];

        Kernel_Elec::ker_from_rho(K1, rho_ele, 1, 0, Dimension::F);
        Kernel_Elec::ker_from_rho(K2, rho_ele, 1, 0, Dimension::F);
    }
    return 0;
}

};  // namespace PROJECT_NS

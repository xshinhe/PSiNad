#include "Kernel_Elec_MMD.h"

#include "Kernel_Dimension.h"
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

/**
 * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
 */
int Kernel_Elec_MMD::rho_focus(num_complex* rho, int iocc, double gamma_ou, double gamma_uu, int fdim, bool rand_act,
                               bool pure_phase, bool cont_phase) {
    double randu;
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
    switch (mmd_type) {
        case MMDPolicy::MMF:
            gamma_uu = (std::sqrt(scale * Kernel_Dimension::F + 1) - 1) / Kernel_Dimension::F;
            gamma_ou = sqrt(gamma_uu * (1.0f + gamma_uu));
            break;
        case MMDPolicy::TWA:
            gamma_uu = 0.0f;
            gamma_ou = sqrt(0.5f * scale);
            break;
        case MMDPolicy::MID: {
            double scale_k = PM->get<num_real>("scale_k", LOC(), 1.0);
            double r       = (Kernel_Dimension::F + 2 * scale_k) / (1 + scale_k);
            gamma_uu       = (std::sqrt(r / (1 + scale_k) * scale + 1) - 1) / r;
            gamma_ou       = std::sqrt((1 + scale_k) * gamma_uu * (1 + gamma_uu));
            break;
        }
    }
    double gamma_uu_in = PM->get<num_real>("gamma_uu", LOC(), -1.0);
    double gamma_ou_in = PM->get<num_real>("gamma_ou", LOC(), -1.0);
    if (gamma_uu_in > 0 && gamma_ou_in > 0) {
        gamma_uu = gamma_uu_in;
        gamma_ou = gamma_ou_in;
    }
    pure_phase = PM->get<bool>("pure_phase", LOC(), true);
    cont_phase = PM->get<bool>("cont_phase", LOC(), true);
    rand_act   = PM->get<bool>("rand_act", LOC(), false);
}

void Kernel_Elec_MMD::init_calc_impl(int stat) {
    Kernel_Elec::w[0] = (rand_act) ? num_complex(Kernel_Dimension::F) : 1.0e0;

    *Kernel_Elec::occ_ele = Kernel_Elec::occ0;                                                   // useless
    *Kernel_Elec::occ_nuc = Kernel_Elec::occ0;                                                   // useless
    rho_focus(Kernel_Elec::rho_ele, Kernel_Elec::occ0, gamma_ou, gamma_uu, Kernel_Dimension::F,  //
              rand_act, pure_phase, cont_phase);
    Kernel_Elec::ker_from_rho(Kernel_Elec::rho_nuc, Kernel_Elec::rho_ele, 1, 0, Kernel_Dimension::F);

    Kernel_Elec::ker_from_rho(Kernel_Elec::K0, Kernel_Elec::rho_ele, 1, 0, Kernel_Dimension::F);
}

int Kernel_Elec_MMD::exec_kernel_impl(int stat) {
    Kernel_Elec::ker_from_rho(Kernel_Elec::Kt, Kernel_Elec::rho_ele, 1, 0, Kernel_Dimension::F);
    return 0;
}

};  // namespace PROJECT_NS

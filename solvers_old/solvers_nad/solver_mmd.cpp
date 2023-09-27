#include "solver_mmd.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

MMD_Solver::MMD_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    std::string phase_flag = Param_GetT(std::string, parm, "phase_flag", "#duni");
    phase_type             = samp_phase::_dict.at(phase_flag);
    std::string focus_flag = Param_GetT(std::string, parm, "focus_flag", "#mmf");
    focus_type             = focus_action::_dict.at(focus_flag);

    Param_GetV(Fref, iparm, F);
    Param_GetV(scale, iparm, 1.0f);
    Param_GetV(sums, iparm, 0);

    if (sums == 1) {
        tcf_weight = phys::math::iu * num_real(F);
        rand_catalog(&occrand, 1, true, 0, F - 1);
    }  // tcf_weight = 1 instead

    switch (focus_type) {  // gamma_uu is nothing but gamma0
        case focus_action::mmf:
            gamma_uu = GAMMA_WIGNER(Fref);
            gamma_uu *= scale;
            gamma_ou = sqrt(gamma_uu * (1.0f + gamma_uu));
            break;
        case focus_action::twa:
            gamma_uu = 0.0f;
            gamma_ou = sqrt(0.5f * (1 - (Fref - 2) * gamma_uu * gamma_uu));
            gamma_uu *= scale;
            gamma_ou *= sqrt(scale);
            break;
        case focus_action::ddd:
            Param_GetV(gamma_uu, iparm, GAMMA_WIGNER(Fref));
            Param_GetV(gamma_ou, iparm, sqrt(0.5f * (1 - (Fref - 2) * gamma_uu * gamma_uu)));  // @bugs
            gamma_uu *= scale;
            gamma_ou *= scale;
            LOG(FATAL);
            break;
        case focus_action::mms1: {
            num_real randu;
            rand_uniform(&randu);
            gamma_uu = randu * (-1.5f + sqrt(2.25f + 3 * Fref)) / Fref;
            gamma_ou = sqrt(gamma_uu * (1 + gamma_uu));
            break;
        }
        case focus_action::mms2: {
            num_real randu, randu2;
            rand_uniform(&randu), rand_uniform(&randu2);
            if (randu2 > randu) randu = randu2;
            gamma_uu = randu * (-4.0f / 3.0f + sqrt(16.0f / 9.0f + 2 * Fref)) / Fref;
            gamma_ou = sqrt(gamma_uu * (1 + gamma_uu));
            break;
        }
        case focus_action::mms3: {
            num_real randu;
            rand_uniform(&randu);
            gamma_uu = -log(1.0f - randu) / (1.0f + sqrt(1.0f + 2 * Fref));
            gamma_ou = sqrt(gamma_uu * (1 + gamma_uu));
            break;
        }
        default:
            LOG(FATAL);
    }
    eom_type = elec_eom::rho;  // use rho for EOM

    std::string suffix;
    suffix += phase_flag;
    suffix += focus_flag;
    suffix += std::to_string(scale);
    save = name() + suffix + "_" + pForceField->tag;
};

MMD_Solver::~MMD_Solver(){};

int MMD_Solver::init(const int& itraj) {
    CHECK_EQ(eom_type, elec_eom::rho);
    NadTraj_Solver::init(itraj);

    if (ini_type == elec_init::occ || ini_type == elec_init::rot || ini_type == elec_init::def ||
        ini_type == elec_init::eig || ini_type == elec_init::d2a) {  // ignore old values and resample in this case
        num_real randu;
        switch (phase_type) {
            case samp_phase::duni:
                for (int j = 1; j < F; ++j) {
                    rand_uniform(&randu);
                    randu          = phys::math::twopi * randu;
                    rho[0 * F + j] = cos(randu) + phys::math::im * sin(randu);
                    rho[j * F + 0] = CONJ_OF(rho[0 * F + j]);
                }
                for (int i = 1; i < F; ++i) {
                    for (int j = i + 1; j < F; ++j) {
                        rho[i * F + j] = rho[0 * F + j] / rho[0 * F + i];
                        rho[j * F + i] = CONJ_OF(rho[i * F + j]);
                    }
                }
                break;
            case samp_phase::funi:
                for (int i = 0; i < F; ++i) {
                    for (int j = i + 1; j < F; ++j) {
                        rand_uniform(&randu);
                        randu          = phys::math::twopi * randu;
                        rho[i * F + j] = cos(randu) + phys::math::im * sin(randu);
                        rho[j * F + i] = CONJ_OF(rho[i * F + j]);
                    }
                }
                break;
            case samp_phase::ddis:
                for (int j = 1; j < F; ++j) {
                    rand_uniform(&randu);
                    randu          = phys::math::halfpi * (int(randu / 0.25f) + 1);
                    rho[0 * F + j] = cos(randu) + phys::math::im * sin(randu);
                    rho[j * F + 0] = CONJ_OF(rho[0 * F + j]);
                }
                for (int i = 1; i < F; ++i) {
                    for (int j = i + 1; j < F; ++j) {
                        rho[i * F + j] = rho[0 * F + j] / rho[0 * F + i];
                        rho[j * F + i] = CONJ_OF(rho[i * F + j]);
                    }
                }
                break;
            case samp_phase::fdis:
                for (int i = 0; i < F; ++i) {
                    for (int j = i + 1; j < F; ++j) {
                        rand_uniform(&randu);
                        randu          = phys::math::halfpi * (int(randu / 0.25f) + 1);
                        rho[i * F + j] = cos(randu) + phys::math::im * sin(randu);
                        rho[j * F + i] = CONJ_OF(rho[i * F + j]);
                    }
                }
                break;
            default:
                LOG(FATAL);
        }

        if (sums == 1) {
            rand_catalog(&occrand);
        } else {
            occrand = occ0;
        }

        for (int i = 0; i < F; ++i) {
            for (int j = 0; j < F; ++j) {
                if (i == j) {
                    rho[i * F + j] = (i == occrand) ? phys::math::iu : phys::math::iz;
                } else if (i == occrand || j == occrand) {
                    rho[i * F + j] *= gamma_ou;
                } else {
                    rho[i * F + j] *= gamma_uu;
                }
            }
        }
        if (ini_type == elec_init::d2a) {
            ARRAY_MATMUL_TRANS1(workc, T, rho, F, F, F);
            ARRAY_MATMUL(rho, workc, T, F, F, F);
        }
        for (int i = 0; i < FF; ++i) gmat[i] = phys::math::iz;
    }
    return 0;
}
int MMD_Solver::kernel0(num_complex* rhox, const int& F) {
    for (int i = 0; i < FF; ++i) rhox[i] = rho[i];
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    // ARRAY_SHOW(rhox, F, F);
    return 0;
}
int MMD_Solver::kernelt(num_complex* rhox, const int& F) { return kernel0(rhox, F); }

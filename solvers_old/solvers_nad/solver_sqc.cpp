#include "solver_sqc.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

/**
 * @brief window sampling and binning in sqc
 */
namespace window {
int samp_mvc_window(num_real* mvc, int& iocc, const int& fdim, const int& window_type) {
    const num_real gm0 = GAMMA_WIGNER(2);
    switch (window_type) {
        case window::tri: {
            num_real tmp2[2];
            rand_uniform(tmp2, 2);
            while (tmp2[0] + tmp2[1] > 1.0f) rand_uniform(tmp2, 2);
            mvc[iocc] = tmp2[0];
            tmp2[1]   = 1.0f - mvc[iocc];
            for (int i = 0; i < fdim; ++i) {
                if (i != iocc) {
                    rand_uniform(tmp2, 1);
                    mvc[i] = tmp2[0] * tmp2[1];
                }
            }
            mvc[iocc] += 1;
            break;
        }
        case window::sqr: {
            rand_uniform(mvc, fdim, 2.0f * gm0);
            mvc[iocc] += 1;
            break;
        }
        default:
            LOG(FATAL);
    }
    samp_mvc_focus(mvc, fdim);
    return 0;
}

int rho_eac(num_complex* rho, num_complex* eac, const int& fdim, const int& window_type) {
    const num_real gm0 = GAMMA_WIGNER(2.0f), gm1 = 1 + gm0, gmh = 0.5f + gm0;

    // normalize rhonm=1 with phase factor
    ARRAY_OUTER_TRANS2(rho, eac, eac, fdim, fdim);
    for (int i = 0; i < fdim * fdim; ++i) rho[i] = rho[i] / ABS_OF(rho[i]);

    // then set zeros (quantize to zero by window function)
    switch (window_type) {
        case window::tri: {
            for (int i = 0; i < fdim; ++i) {
                int idx1 = i * (fdim + 1);
                int idx2 = idx1;  // for (i,j) & (j,i)
                for (int k = 0; k < fdim; ++k) {
                    if ((k != i && NORM_OF(eac[k]) > 1) || (k == i && NORM_OF(eac[k]) < 1)) {
                        rho[idx1] = phys::math::iz;
                        break;
                    }
                }
                idx1++;        // (i,j++)
                idx2 += fdim;  // (i++,j)
                for (int j = i + 1; j < fdim; ++j) {
                    for (int k = 0; k < fdim; ++k) {
                        if ((k != i && NORM_OF(eac[k]) > 1) || (k == i && NORM_OF(eac[k]) < 0.5f) ||
                            (k == j && NORM_OF(eac[k]) < 0.5f)) {
                            rho[idx1] = phys::math::iz, rho[idx2] = phys::math::iz;
                            break;
                        }
                    }
                    idx1++;        // (i,j++)
                    idx2 += fdim;  // (i++,j)
                }
            }
            break;
        }
        case window::sqr: {
            for (int i = 0; i < fdim; ++i) {
                int idx1 = i * (fdim + 1);
                int idx2 = idx1;  // for rho(j,i)
                for (int k = 0; k < fdim; ++k) {
                    if ((k != i && std::abs(NORM_OF(eac[k]) - gm0) > gm0) ||
                        (k == i && std::abs(NORM_OF(eac[k]) - gm1) > gm0)) {
                        rho[idx1] = phys::math::iz;
                        break;
                    }
                }
                idx1++;        // (i,j++)
                idx2 += fdim;  // (i++,j)
                for (int j = i + 1; j < fdim; ++j) {
                    for (int k = 0; k < fdim; ++k) {
                        if ((k != i && std::abs(NORM_OF(eac[k]) - gm0) > gm0) ||
                            (k == i && std::abs(NORM_OF(eac[k]) - gmh) > gm0) ||
                            (k == j && std::abs(NORM_OF(eac[k]) - gmh) > gm0)) {
                            rho[idx1] = phys::math::iz, rho[idx2] = phys::math::iz;
                            break;
                        }
                    }
                    idx1++;        // (i,j++)
                    idx2 += fdim;  // (i++,j)
                }
            }
            break;
        }
        default:
            LOG(FATAL);
    }
    return 0;
}
};  // namespace window



SQC_Solver::SQC_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    window_type = window::_dict.at(Param_GetT(std::string, parm, "window_type", "#tri"));

    switch (window_type) {
        case window::tri:
            gamma0 = 1.0f / 3.0f;
            gammat = gamma0;
            break;
        case window::sqr:
            gamma0 = 0.5f;
            gammat = gamma0;
            break;
        default:
            LOG(FATAL);
    }
    save = name() + "_" + pForceField->tag;
};

SQC_Solver::~SQC_Solver(){};

int SQC_Solver::init_occ2eac(const int& itraj) {
    window::samp_mvc_window(mvc, occ0, F, window_type);
    eac_mvc(eac0, mvc, F);
    return 0;
}
int SQC_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);
    if (eom_type != elec_eom::eaccv) {  // override gmat
        for (int i = 0, idx = 0; i < F; ++i)
            for (int j = 0; j < F; ++j, ++idx) gmat[idx] = (i == j) ? gamma0 : phys::math::iz;
    }
    return 0;
}
int SQC_Solver::kernel0(num_complex* rhox, const int& F) { return window::rho_eac(rhox, eac, F, window_type); }
int SQC_Solver::kernelt(num_complex* rhox, const int& F) { return window::rho_eac(rhox, eac, F, window_type); }

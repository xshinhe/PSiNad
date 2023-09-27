#include "solver_cmm.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

CMM_Solver::CMM_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    xi0 = 1.0f;
    Param_GetV(gamma0, iparm, GAMMA_WIGNER(F));
    gammat     = (1 - gamma0) / (1.0f + F * gamma0);
    totact     = (1 + F * gamma0) / xi0;
    xit        = (1 + F * gammat) / totact;
    tcf_weight = num_complex(F);

    std::string suffix;
    suffix += ((rep_type == representation::adiabatic) ? "adia" : "");
    suffix += ((eom_type == elec_eom::eaccv) ? "cv" : "");
    suffix += std::to_string(gamma0);

    save = name() + suffix + "_" + pForceField->tag;
};

CMM_Solver::~CMM_Solver(){};

int CMM_Solver::init_occ2eac(const int& itraj) {
    samp_mvc_sphere(mvc, 2.0f * totact, F);
    eac_mvc(eac0, mvc, F);
    return 0;
}

int CMM_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);
    if (eom_type != elec_eom::eaccv) {  // override gmat
        for (int i = 0, idx = 0; i < F; ++i)
            for (int j = 0; j < F; ++j, ++idx) gmat[idx] = (i == j) ? gamma0 : phys::math::iz;
    }
    return 0;
}
int CMM_Solver::kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F) {
    rho_eac(rhox, eac, F);
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    for (int idx = 0, i = 0; i < F; ++i)
        for (int j = 0; j < F; ++j, ++idx)
            rhox[idx] = xic * rhox[idx] - ((i == j) ? gammac * phys::math::iu : phys::math::iz);
    return 0;
}

int CMM_Solver::kernel0(num_complex* rhox, const int& F) { return kernel_cmm(rhox, xi0, gamma0, F); }
int CMM_Solver::kernelt(num_complex* rhox, const int& F) { return kernel_cmm(rhox, xit, gammat, F); }

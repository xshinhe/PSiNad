#include "solver_lsc.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

LSC_Solver::LSC_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    gamma0             = Param_GetV(gamma0, iparm, 0.250f);  // default 0.25
    gammat             = Param_GetV(gammat, iparm, gamma0);  // default 0.25
    sigma0             = std::sqrt(gamma0);
    sigmat             = std::sqrt(gammat);
    std::string gsflag = Param_GetT(std::string, iparm, "gsflag", "0");  // where to sample gauss, "0"/"t"
    variance           = sigma0 * sigma0;
    tcf_weight         = phys::math::iu / (variance * variance);

    if (gsflag == "0") {
        variance = sigma0 * sigma0;
    } else if (gsflag == "t") {
        variance = sigmat * sigmat;
    }

    std::string suffix = std::to_string(gamma0);

    save = name() + suffix + "_" + pForceField->tag;
};

LSC_Solver::~LSC_Solver(){};

int LSC_Solver::init_occ2eac(const int& itraj) {
    samp_mvc_gauss(mvc, variance, F);
    eac_mvc(eac0, mvc, F);
    return 0;
}

int LSC_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);
    if (eom_type != elec_eom::eaccv) {  // override gmat
        for (int i = 0, idx = 0; i < F; ++i)
            for (int j = 0; j < F; ++j, ++idx) gmat[idx] = (i == j) ? gamma0 : phys::math::iz;  // gamma0 = 0.5?
    }
    return 0;
}

int LSC_Solver::kernel_lsc(num_complex* rhox, num_real& gammac, const int& F) {
    rho_eac(rhox, eac, F);
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    for (int idx = 0, i = 0; i < F; ++i)
        for (int j = 0; j < F; ++j, ++idx)
            rhox[idx] = rhox[idx] - ((i == j) ? gammac * phys::math::iu : phys::math::iz);
    return 0;
}

int LSC_Solver::kernel0(num_complex* rhox, const int& F) { return kernel_lsc(rhox, gamma0, F); }
int LSC_Solver::kernelt(num_complex* rhox, const int& F) { return kernel_lsc(rhox, gammat, F); }

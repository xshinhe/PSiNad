#include "solver_smm.h"

#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

SMM_Solver::SMM_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
    ALLOCATE_PTR_TO_VECTOR(xi1, FF);

    tcf_weight = num_complex(F);
    // tcf_weight = num_complex(1.0f);

    Param_Reset(eom_type, elec_eom::rho);

    save = name() + "_" + pForceField->tag;
};

SMM_Solver::~SMM_Solver(){};

int SMM_Solver::init_occ2eac(const int& itraj) {
    rand_sphere(xi1, FF - 1, sqrt((FF - 1.0f) / F));
    xi1[FF - 1] = sqrt(1.0f / F);

    int Fadd1      = F + 1;
    num_real save1 = 0.0f, save2 = 0.0f;

    for (int i = 1; i < F; ++i) xi1[(i - 1) * Fadd1] /= sqrt(i * i + i);
    xi1[FF - 1] /= sqrt(1.0f * F);

    for (int i = 0; i < F; ++i) {
        save1 = xi1[i * Fadd1];
        xi1[i * Fadd1] -= i * save2;
        save2 = save1;
        for (int j = i + 1; j < F; ++j) xi1[i * Fadd1] += xi1[j * Fadd1];
    }
    for (int i = 0, idx = 0; i < F; ++i) {
        rho0[i * Fadd1] = xi1[i * Fadd1];
        for (int j = i + 1; j < F; ++j, ++idx) {
            rho0[i * F + j] = phys::math::sqrthalf * (xi1[i * F + j] + phys::math::im * xi1[j * F + i]);
            rho0[j * F + i] = phys::math::sqrthalf * (xi1[i * F + j] - phys::math::im * xi1[j * F + i]);
        }
    }

    // focus operation
    // for (int i = 0, idxd = 0; i < F; ++i, idxd += Fadd1) rho0[idxd] = (i == occ0) ? phys::math::iu : phys::math::iz;

    for (int i = 0; i < FF; ++i) rho[i] = rho0[i];

    return 0;
}

int SMM_Solver::init(const int& itraj) {
    NadTraj_Solver::init(itraj);
    return 0;
}
int SMM_Solver::kernel_scmm(num_complex* rhox, const int& F) {
    for (int i = 0; i < FF; ++i) rhox[i] = rho[i];
    if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
        ARRAY_MATMUL(workc, T, rhox, F, F, F);
        ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
    }
    return 0;
}
int SMM_Solver::kernel0(num_complex* rhox, const int& F) { return kernel_scmm(rhox, F); }
int SMM_Solver::kernelt(num_complex* rhox, const int& F) { return kernel_scmm(rhox, F); }

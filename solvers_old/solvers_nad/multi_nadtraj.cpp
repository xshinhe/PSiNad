#include "multi_nadtraj.h"

#include "../utils/definitions.h"
#include "../utils/nad_utils.h"

using namespace ARRAY_EG;

Multi_NadTraj_Solver::Multi_NadTraj_Solver(const Param& iparm, Model* pM, const int& _M) : NadTraj_Solver(iparm, pM) {
    CHECK_GT(_M, 0);
    M  = _M;
    MN = M * N, MF = M * F, MNN = M * NN, MFF = M * FF, MNFF = M * NFF;

    try {
        ALLOCATE_PTR_TO_VECTOR(nrs, MN);
        ALLOCATE_PTR_TO_VECTOR(nps, MN);
        ALLOCATE_PTR_TO_VECTOR(nms, MN);
        ALLOCATE_PTR_TO_VECTOR(nfs, MN);
        ALLOCATE_PTR_TO_VECTOR(vpeses, M);
        ALLOCATE_PTR_TO_VECTOR(grads, MN);
        // ALLOCATE_PTR_TO_VECTOR(hesses, MNN);
        ALLOCATE_PTR_TO_VECTOR(Vs, MFF);
        ALLOCATE_PTR_TO_VECTOR(dVs, MNFF);
        // ALLOCATE_PTR_TO_VECTOR(ddVs, M*NNFF);

        ALLOCATE_PTR_TO_VECTOR(Es, MF);
        ALLOCATE_PTR_TO_VECTOR(Ts, MFF);
        ALLOCATE_PTR_TO_VECTOR(dEs, MNFF);
        // ALLOCATE_PTR_TO_VECTOR(nacvs, MNFF);
        ALLOCATE_PTR_TO_VECTOR(Ls, MF);
        ALLOCATE_PTR_TO_VECTOR(Ss, MFF);
        ALLOCATE_PTR_TO_VECTOR(rhos, MFF);
        ALLOCATE_PTR_TO_VECTOR(Us, MFF);
    } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }
}
Multi_NadTraj_Solver::~Multi_NadTraj_Solver(){};

int Multi_NadTraj_Solver::ff_calc1(const int& level, const bool& refered) { return multi_ff_calc1(level, refered); }
int Multi_NadTraj_Solver::multi_ff_calc1(const int& level,
                                         const bool& refered) {  // multi-forcefield calculation at a fixed nr, np
    int succ = 0;
    for (int mu = 0; mu < M; ++mu) {
        num_real* vpes_mu = vpeses + mu;
        num_real* grad_mu = grads + mu * N;
        num_real* nr_mu   = nrs + mu * N;
        num_real* np_mu   = nps + mu * N;
        num_real* nm_mu   = nms + mu * N;
        num_real* E_mu    = Es + mu * F;  ///< solve diabatic propagators
        num_real* T_mu    = Ts + mu * FF;
        num_real* L_mu    = Ls + mu * F;  ///< solve adiabatic propagators
        num_complex* S_mu = Ss + mu * FF;
        num_complex* U_mu = Us + mu * FF;
        num_real* V_mu    = Vs + mu * FF;
        num_real* dV_mu   = dVs + mu * NFF;
        num_real* dE_mu   = dEs + mu * NFF;

        switch (pForceField->type) {
            case ForceFieldModel: {
                if (succ == 0)
                    succ = pForceField->ForceField_npes(vpes_mu, grad_mu, hess, nr_mu, np_mu, level, N, itraj, istep);
                if (succ == 0) succ = pForceField->ForceField_epes(V_mu, dV_mu, ddV, nr_mu, level, N, F, itraj, istep);
                break;
            }
            case ForceFieldOnTheFly: {
                if (rep_type != representation::onthefly) Param_Reset(rep_type, representation::onthefly);
                if (succ == 0) succ = pForceField->ForceField_epes(E_mu, dE_mu, ddE, nr_mu, level, N, F, itraj, istep);
                break;
            }
            default:
                LOG(FATAL);
        }
        solve_transform(H, dH, ddH, S_mu, L_mu, dL, ddL, T_mu, E_mu, dE_mu, ddE, V_mu, dV_mu, ddV, nr_mu, np_mu, nm_mu,
                        rep_type, level, N, F, workr, workc, refered);
        solve_Ut(U_mu, S_mu, L_mu, T_mu, E_mu, dt, rep_type, N, F, workr, workc);
    }
    return succ;
}
int Multi_NadTraj_Solver::ff_calc2() { return multi_ff_calc2(); }
int Multi_NadTraj_Solver::multi_ff_calc2() {
    LOG(FATAL);
    return 0;
}
int Multi_NadTraj_Solver::evolve_elec(num_complex* Uevolve) { return multi_evolve_elec(Uevolve); }
int Multi_NadTraj_Solver::multi_evolve_elec(num_complex* Uevolve) {
    LOG(FATAL);
    return 0;
}
int Multi_NadTraj_Solver::update_p(const num_real& dt_in) {
    for (int mu = 0, idx = 0; mu < M; ++mu)
        for (int i = 0; i < N; ++i, ++idx) nps[idx] -= nfs[idx] * dt_in;
    return 0;
}
int Multi_NadTraj_Solver::update_r(const num_real& dt_in) {
    for (int mu = 0, idx = 0; mu < M; ++mu)
        for (int i = 0; i < N; ++i, ++idx) nrs[idx] += nps[idx] / nms[idx] * dt_in;
    return 0;
}
int Multi_NadTraj_Solver::update_thermo(const num_real& dt_in) {
    LOG(FATAL);  // not implement yet
    return 0;
}
int Multi_NadTraj_Solver::init(const int& itraj) { return multi_init(itraj); };
int Multi_NadTraj_Solver::multi_init(const int& itraj) {
    LOG(FATAL);
    return 0;
};

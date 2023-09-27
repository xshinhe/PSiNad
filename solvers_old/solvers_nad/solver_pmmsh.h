#ifndef PMMSH_SOLVER_H
#define PMMSH_SOLVER_H

class PMMSH_Solver : public PMM_Solver {
    PMMSH_Solver(const Param& iparm, Model* pM) : PMM_Solver(iparm, pM){};
    virtual int multi_ff_calc2() {
        num_real* ForceMats  = (rep_type == representation::diabatic) ? dVs : dEs;
        num_real* ForceMat_0 = ForceMats;
        num_real* nf_0       = nfs;
        num_real* grad_0     = grads;
        num_complex* rho_0   = rhos;

        pForceField->reduce_force(fmean, rho_0, ForceMat_0, N, F);
        for (int i = 0; i < N; ++i) nf_0[i] = grad_0[i] + fmean[i];  // total force

        // calculate correllated force
        for (int mu = 1; mu < M; ++mu) {
            num_real* ForceMat_mu   = ForceMats + mu * NFF;
            num_real* nf_mu         = nfs + mu * N;
            num_real* grad_mu       = grads + mu * N;
            num_complex* diffrho_mu = diffrhos + mu * FF;
            for (int i = 0; i < N; ++i) nf_mu[i] = grad_mu[i];
            pForceField->reduce_force(fmean, rho_0, ForceMat_mu, N, F);
            for (int i = 0; i < N; ++i) nf_mu[i] += fmean[i];
            pForceField->reduce_force(fmean, diffrho_mu, ForceMat_0, N, F);
            for (int i = 0; i < N; ++i) nf_mu[i] += fmean[i];
        }
        return 0;
    }
    virtual ~PMMSH_Solver();

    static inline std::string name() { return "pmmsh"; }
};

#endif  // PMMSH_SOLVER_H
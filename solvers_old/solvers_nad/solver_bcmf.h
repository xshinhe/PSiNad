#ifndef BCMF_SOLVER_H
#define BCMF_SOLVER_H
#include "nadtraj.h"

// native adiabatic ???

class BCMF_Solver : public NadTraj_Solver {
   public:
    BCMF_Solver(Param iparm, Model* pM) : NadTraj_Solver(iparm, pM) {
        save = BCMF_Solver::name() + "_" + pForceField->tag;
        Param_Reset(rep_type, representation::adiabatic);
        Param_Reset(eom_type, elec_eom::eac);

        Param_AutoDefault(int, bc_flag, 1, iparm);  // do bifurcation, otherwise only phase correction

        try {
            ALLOCATE_PTR_TO_VECTOR(scale_arr, F);
            ALLOCATE_PTR_TO_VECTOR(reflect_arr, F);
        } catch (std::bad_alloc& e) { LOG(FATAL) << e.what(); }
    };

    virtual ~BCMF_Solver(){};

    static inline std::string name() { return "bcmf"; }

    virtual int init_occ2eac(const int& itraj) { return NadTraj_Solver::init_occ2eac(itraj); }

    virtual int init(const int& itraj) {
        NadTraj_Solver::init(itraj);  // solver eigen problem here
        for (int i = 0; i < FF; ++i) gmat[i] = phys::math::iz;

        if (ini_type == elec_init::d2a) {  //? @TODO
            double rand_tmp, sum = 0.0f;
            rand_uniform(&rand_tmp);
            for (int i = 0; i < F; ++i) {
                sum += NORM_OF(eac[i]);
                if (rand_tmp < sum) {
                    occt = i;  // project from diabatic occ0 to adiabatic occt
                    break;
                }
            }
        }
        return 0;
    }
    virtual int kernel0(num_complex* rhox, const int& F) {
        rho_eac(rhox, eac, F);
        if (rep_type == representation::adiabatic && ini_type == elec_init::d2a) {
            ARRAY_MATMUL(workc, T, rhox, F, F, F);
            ARRAY_MATMUL_TRANS2(rhox, workc, T, F, F, F);
        }
        return 0;
    }
    virtual int kernelt(num_complex* rhox, const int& F) { return kernel0(rhox, F); }

    virtual int ff_calc1(const int& level = 1) {  // multi-forcefield calculation at a fixed nr, np
        int succ = 0;
        // NOTE: H is just adiabtic H, but U is of Heff (with phase corrected)
        rho_eac(rho, eac, F);
        switch (pForceField->type) {
            case ForceFieldModel:
                pForceField->ForceField_npes(vpes, grad, hess, nr, np, level, N);
                succ = pForceField->ForceField_epes(V, dV, ddV, nr, level, N, F);
                if (succ == 0) {
                    solve_transform_correctphase(H, rho, S, L, dL, ddL, T, E, dE, ddE, V, dV, ddV, nr, np, nm, rep_type,
                                                 level, N, F, workr, workc);
                } else
                    LOG(WARNING) << "failed to call forcefield";
                break;
            case ForceFieldOnTheFly:
                rep_type = -1;  // rep_type = -1 for on-the-fly adiabatic
                succ     = pForceField->ForceField_epes(E, dE, ddE, nr, level, N, F);
                if (succ == 0) {
                    solve_transform_correctphase(H, rho, S, L, dL, ddL, T, E, dE, ddE, V, dV, ddV, nr, np, nm, rep_type,
                                                 level, N, F, workr, workc);
                } else
                    LOG(WARNING) << "failed to call forcefield";
                break;
            default:
                LOG(FATAL) << "Unknown ForceField type";
                break;
        }

        if (bc_flag < 1) return succ;

        // otherwise do bifurcation

        succ = ff_calc2();

        // mean WP (effective WP)
        rho_eac(rho, eac, F);
        double Eavg = REAL_OF(ARRAY_TRACE2(rho, H, F, F));  // @CHECK???
        double Ekin = 0.0f, fdotp = 0.0f, fnorm = 0.0f;
        for (int j = 0; j < N; ++j) {
            Ekin += 0.5f * np[j] * np[j] / nm[j];
            fdotp += nf[j] * np[j];
            fnorm += nf[j] * nf[j];
        }
        bool reflect_mean = ((fdotp / fnorm) < dt);

        // WPS
        for (int i = 0; i < F; ++i) {
            // determine scale factor
            double s2    = 1.0 + (Eavg - REAL_OF(H[i * (F + 1)])) / Ekin;
            scale_arr[i] = (s2 > 0) ? sqrt(s2) : 10e-12;  // unaccessiable
            // if reflected
            fdotp = 0.0f, fnorm = 0.0f;
            for (int j = 0; j < N; ++j) {
                double nfj = grad[j] + dE[j * FF + i * (F + 1)] - fcorr[j];
                fdotp += nfj * np[j] * scale_arr[i];
                fnorm += nfj * nfj;
            }
            reflect_arr[i] = ((fdotp / fnorm) < dt);
        }

        // divide into RG, NRG
        double refpop = 0.0f;
        for (int i = 0; i < F; ++i)
            if (reflect_arr[i]) refpop += NORM_OF(eac[i]);

        // first and second kind
        if (0.05 < refpop && refpop < 0.95) {  // bifurcation of WPs
            double rand_tmp;
            rand_uniform(&rand_tmp);
            LOG(WARNING) << "jump with rand=" << rand_tmp << ", in (" << refpop << "," << 1 - refpop << ")";
            if (rand_tmp < refpop) {
                for (int i = 0; i < F; ++i) {
                    if (reflect_arr[i])
                        eac[i] /= sqrt(refpop);
                    else
                        eac[i] = phys::math::iz;
                }
            } else {
                for (int i = 0; i < F; ++i) {
                    if (!reflect_arr[i]) {
                        eac[i] /= sqrt(1.0f - refpop);
                    } else {
                        eac[i] = phys::math::iz;
                    }
                }
            }
            rho_eac(rho, eac, F);
            double Eavgp = REAL_OF(ARRAY_TRACE2(rho, H, F, F));  // @CHECK???
            double s2    = 1.0 + (Eavg - Eavgp) / Ekin;
            double s     = (s2 > 0) ? sqrt(s2) : 10e-12;  // unaccessiable ??
            for (int j = 0; j < N; ++j) np[j] *= s;

            // exception??????
        } else if (reflect_mean && refpop >= 1) {  // remove unaccessiable state
            LOG(WARNING) << "jump third type";

            double norm = 0.0f;
            for (int i = 0; i < F; ++i) {
                if (scale_arr[i] < 10E-8)
                    eac[i] = phys::math::iz;
                else
                    norm += NORM_OF(eac[i]);
            }
            for (int i = 0; i < F; ++i) eac[i] /= sqrt(norm);
            rho_eac(rho, eac, F);
            double Eavgp = REAL_OF(ARRAY_TRACE2(rho, H, F, F));  // @CHECK???
            double s2    = 1.0 + (Eavg - Eavgp) / Ekin;
            double s     = (s2 > 0) ? sqrt(s2) : 10e-12;  // unaccessiable ??
            for (int j = 0; j < N; ++j) np[j] *= s;
        } else {
            // do nothing
        }
        return succ;
    }

    DEFINE_POINTER(num_real, scale_arr);
    DEFINE_POINTER(bool, reflect_arr);
    int bc_flag = 0;
};

#endif  // BCMF_SOLVER_H

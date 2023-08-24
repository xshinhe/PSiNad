#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_GWP final : public Kernel {
   public:
    static int Snuc(num_complex* S,   // [P,P]
                    num_real* x,      // [P,N]
                    num_real* p,      // [P,N]
                    num_real* alpha,  // [N]
                    int P, int N) {
        for (int m = 0, mj = 0, mi = 0, mn = 0; m < P; ++m) {
            for (int n = 0, nj = 0, ni = 0; n < P; ++n, ++mn) {
                num_complex term = 0.0e0;
                for (int j = 0; j < N; ++j, ++mj, ++nj) {
                    term += -0.25 * alpha[j] * (x[mj] - x[nj]) * (x[mj] - x[nj])   //
                            - 0.25 / alpha[j] * (p[mj] - p[nj]) * (p[mj] - p[nj])  //
                            - 0.5 * (p[mj] + p[nj]) * (x[mj] - x[nj]);
                }
                S[mn] = std::exp(term);
            }
        }
        return 0;
    }

    static int Sele(num_complex* S,  // [P,P]
                    num_complex* c,  // [P,F]
                    int P, int F) {
        for (int m = 0, mj = 0, mi = 0, mn = 0; m < P; ++m) {
            for (int n = 0, nj = 0, ni = 0; n < P; ++n, ++mn) {
                num_complex term = 0.0e0;
                for (int i = 0; i < F; ++i, ++mi, ++ni) { term += std::conj(c[mi]) * c[ni]; }
                S[mn] = term;
            }
        }
        return 0;
    }

    static int dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                        num_real* x,            // [P,N]
                        num_real* p,            // [P,N]
                        num_real* m,            // [P,N]
                        num_real* f,            // [P,N]
                        num_real* alpha,        // [N]
                        int P, int N) {
        for (int m = 0, mj = 0, mn = 0; m < P; ++m) {
            for (int n = 0, nj = 0; n < P; ++n, ++mn) {
                num_complex term = 0.0e0;
                for (int j = 0; j < N; ++j, ++mj, ++nj) {
                    term += p[nj] / m[nj] * 0.5 * alpha[j] *
                            (x[mj] - x[nj] + (p[mj] + p[nj]) / (phys::math::im * alpha[j]));
                    term -= f[nj] * 0.5 * ((p[mj] - p[nj]) + phys::math::im * alpha[j] * (x[mj] - x[nj])) / alpha[j];
                    term += phys::math::im * p[nj] * p[nj] / (2 * m[nj]);
                }
                dtlnSnuc[mn] = term;
            }
        }
        return 0;
    }

    static int dtSele(num_complex* dtSele,  // [P,P]
                      num_complex* c,       // [P,F]
                      num_complex* dUdt,    // [P,F,F]
                      int P, int F) {
        for (int m = 0, mj = 0, mn = 0; m < P; ++m) {
            for (int n = 0, nj = 0; n < P; ++n, ++mn) {
                num_complex* cm       = c + m * F;
                num_complex* cn       = c + n * F;
                num_complex* dUdtn    = dUdt + n * F * F;
                num_complex* dtSelemn = dtSele + mn;
                ARRAY_MATMUL3_CONJ1(dtSelemn, cm, dUdt, cn, 1, F, F, 1);
            }
        }
        return 0;
    }

    static int reduce_density(num_complex* rhored,  // [F,F]
                              num_complex* Acoeff,  // [P]
                              num_complex* Snuc,    // [P,P]
                              num_complex* c,       // [P,F]
                              int P, int F) {       //@bug
        Eigen::Map<EigMX<num_complex>> Map_rhored(rhored, Dimension::F, Dimension::F);
        Eigen::Map<EigMX<num_complex>> Map_Acoeff(Acoeff, Dimension::P, Dimension::1);
        Eigen::Map<EigMX<num_complex>> Map_Snuc(Snuc, Dimension::P, Dimension::P);
        Eigen::Map<EigMX<num_complex>> Map_c(c, Dimension::P, Dimension::F);
        Map_rhored = Map_c.adjoint() * Map_Acoeff.asDiagonal().adjoint() * Map_Snuc * Map_Acoeff.asDiagonal() * Map_c;
        return 0;
    }

    static int effective_ham(num_complex* Heff,    // [P,P]
                             num_real* V,          // [P,F,F]
                             num_real* dV,         // [P,N,F,F]
                             num_real* x,          // [P,N]
                             num_real* p,          // [P,N]
                             num_real* m,          // [P,N]
                             num_complex* c,       // [P,F]
                             int P, int N, int F,  // Dimensions
                             int rep_flag) {
        for (int m = 0, mj = 0, mn = 0; m < P; ++m) {
            for (int n = 0, nj = 0; n < P; ++n, ++mn) {
                num_real* Vn    = V + n * FF;
                num_real* Vm    = V + m * FF;
                num_real* dVn   = dV + n * NFF;
                num_real* dVm   = dV + m * NFF;
                num_complex* cn = c + n * F;
                num_complex* cm = c + m * F;

                num_complex Tmn = 0.0e0;
                num_complex Vmn = 0.0e0;

                // ARRAY_MATMUL_CONJ2(rhonm, cn, cm, F, 1, F);
                // Vmn += 0.5e0 * ARRAY_TRACE2(Vn, rhonm, F, F);
                // Vmn += 0.5e0 * ARRAY_TRACE2(Vm, rhonm, F, F);

                for (int i = 0; i < F; ++i) {
                    for (int k = 0; k < F; ++k) Vmn += std::conj(cm[i]) * 0.5 * (Vm[ik] + Vn[ik]) * cn[k];
                }
                for (int j = 0, jik = 0; j < N; ++j, ++mj, ++nj) {
                    num_complex xmnj = 0.5e0 * (x[mj] + x[nj] + (p[mj] - p[nj]) / (phys::math::im * alpha[j]));
                    num_complex pmnj = 0.5e0 * (p[mj] + p[nj] + (phys::math::im * alpha[j]) * (x[mj] - x[nj]));
                    Tmn += alpha[j] / (4 * m[j]) + pmnj * pmnj / (2 * m[j]);
                    for (int i = 0, ik = 0; i < F; ++i) {
                        for (int k = 0; k < F; ++k, ++ik, ++jik) {
                            Vmn += std::conj(cm[i]) * 0.5 * ((xmnj - x[nj]) * dVn[jik] + (xmnj - x[mj]) * dVm[jik]) *
                                   cn[k];
                        }
                    }
                }
                Heff[mn] = Tmn + Vmn;
            }
        }
        return 0;
    };

   private:
    num_complex* Hcoeff;
    num_complex* Acoeff;
    num_complex* c;

    void read_param_impl(Param* PM) {
        gamma1  = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));
        gamma2  = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
        xi1     = (1 + Dimension::F * gamma1);
        xi2     = (1 + Dimension::F * gamma2);
        use_cv  = PM->get<bool>("use_cv", LOC(), false);
        use_wmm = PM->get<bool>("use_wmm", LOC(), false);  // @disable
    }

    void init_data_impl(DataSet* DS) {
        alpha  = DS->reg<num_real>("integrator.alpha", Dimension::N);
        Acoeff = DS->reg<num_complex>("integrator.Acoeff", Dimension::P);
        Hcoeff = DS->reg<num_complex>("integrator.Hcoeff", Dimension::PP);
    }

    void init_calc_impl(int stat = -1) {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            //
        }
        exec_kernel(stat);
    }

    int exec_kernel_impl(int stat) {
        EigenSolve(S, lambdaS, TS, P);
        return 0;
    }
};

};  // namespace PROJECT_NS
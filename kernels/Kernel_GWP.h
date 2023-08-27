#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_GWP final : public Kernel {
   public:
    static int calc_Snuc(num_complex* S,   // [P,P]
                         num_real* x,      // [P,N]
                         num_real* p,      // [P,N]
                         num_real* alpha,  // [N]
                         int P, int N) {
        // Eigen::Map<EigMX<num_complex>> Map_S(S, P, P);
        // Eigen::Map<EigMX<num_real>> Map_x(x, P, N);
        // Eigen::Map<EigMX<num_real>> Map_p(p, P, N);
        // Eigen::Map<EigMX<num_real>> Map_a(alpha, N, 1);

        // auto XX = Map_x * Map_x.transpose();
        // auto XP = Map_x * Map_p.transpose();
        // auto PP = Map_p * Map_p.transpose();

        for (int a = 0, aj = 0, ai = 0, ab = 0; a < P; ++a) {
            for (int b = 0, bj = 0, bi = 0; b < P; ++b, ++ab) {
                num_complex term = 0.0e0;
                for (int j = 0; j < N; ++j, ++aj, ++bj) {
                    term += -0.25 * alpha[j] * (x[aj] - x[bj]) * (x[aj] - x[bj])   //
                            - 0.25 / alpha[j] * (p[aj] - p[bj]) * (p[aj] - p[bj])  //
                            - 0.5 * (p[aj] + p[bj]) * (x[aj] - x[bj]);
                }
                S[ab] = std::exp(term);
            }
        }
        return 0;
    }

    static int calc_Sele(num_complex* S,  // [P,P]
                         num_complex* c,  // [P,F]
                         int P, int F) {
        Eigen::Map<EigMX<num_complex>> Map_S(S, P, P);
        Eigen::Map<EigMX<num_complex>> Map_c(c, P, F);
        Map_S = Map_c.conjugate() * Map_c.transpose();
        return 0;
    }

    static int calc_dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                             num_real* x,            // [P,N]
                             num_real* p,            // [P,N]
                             num_real* m,            // [P,N]
                             num_real* f,            // [P,N]
                             num_real* alpha,        // [N]
                             int P, int N) {
        // Eigen::Map<EigMX<num_complex>> Map_dtlnS(dtlnSnuc, P, P);
        // Eigen::Map<EigMX<num_real>> Map_x(x, P, N);
        // Eigen::Map<EigMX<num_real>> Map_p(p, P, N);
        // Eigen::Map<EigMX<num_real>> Map_m(m, P, N);
        // Eigen::Map<EigMX<num_real>> Map_f(f, P, N);
        // Eigen::Map<EigMX<num_real>> Map_a(alpha, N, 1);
        // auto ve = (Map_p.array() / Map_m.array()).matrix();

        for (int a = 0, aj = 0, ab = 0; a < P; ++a) {
            for (int b = 0, bj = 0; b < P; ++b, ++ab) {
                num_complex term = 0.0e0;
                for (int j = 0; j < N; ++j, ++aj, ++bj) {
                    term += p[bj] / m[bj] * 0.5 * alpha[j] *
                            (x[aj] - x[bj] + (p[aj] + p[bj]) / (phys::math::im * alpha[j]));
                    term += -f[bj] * 0.5 / alpha[j] * ((p[aj] - p[bj]) + phys::math::im * alpha[j] * (x[aj] - x[bj]));
                    term += phys::math::im * p[bj] * p[bj] / (2 * m[bj]);
                }
                dtlnSnuc[ab] = term;
            }
        }
        return 0;
    }

    static int calc_dtlnSele(num_complex* dtlnSele,  // [P,P]
                             num_complex* c,         // [P,F]
                             num_complex* dUdt,      // [P,F,F]
                             int P, int F) {
        // Eigen::Map<EigMX<num_complex>> Map_dtlnSele(dtlnSele, P, P);
        // Eigen::Map<EigMX<num_complex>> Map_c(c, P, F);
        // Eigen::Map<EigMX<num_complex>> Map_c(c, P, F, F);
        // Map_S = Map_c.conjugate() * Map_c.transpose();

        for (int a = 0, aj = 0, ab = 0; a < P; ++a) {
            for (int b = 0, bj = 0; b < P; ++b, ++ab) {
                num_complex* ca    = c + a * F;
                num_complex* cb    = c + b * F;
                num_complex* dUdtb = dUdt + b * F * F;
                dtlnSele[ab]       = ARRAY_INNER_VMV_TRANS1(ca, dUdtb, cb, F, F) / ARRAY_INNER_TRANS1(ca, cb, F);
            }
        }
        return 0;
    }

    static int calc_density(num_complex* rhored,  // [F,F]
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

    static int calc_Hcoeff(num_complex* Hcoeff,  // [P,P]
                           num_real* V,          // [P,F,F] or E [P,F]
                           num_real* dV,         // [P,N,F,F] or dE [P,N,F,F]
                           num_real* x,          // [P,N]
                           num_real* p,          // [P,N]
                           num_real* m,          // [P,N]
                           num_complex* c,       // [P,F]
                           int P, int N, int F,  // Dimensions
                           int rep_flag) {
        int FF  = F * F;
        int NFF = N * FF;

        for (int a = 0, aj = 0, ab = 0; a < P; ++a) {
            for (int b = 0, bj = 0; b < P; ++b, ++ab) {
                num_real* Va    = V + a * FF;
                num_real* Vb    = V + b * FF;
                num_real* dVa   = dV + a * NFF;
                num_real* dVb   = dV + b * NFF;
                num_complex* ca = c + a * F;
                num_complex* cb = c + b * F;

                num_complex Tab = 0.0e0;
                num_complex Vab = 0.0e0;

                Vab += 0.5e0 * ARRAY_INNER_VMV_TRANS1(ca, Va, cb, F, F);
                Vab += 0.5e0 * ARRAY_INNER_VMV_TRANS1(ca, Vb, cb, F, F);

                for (int j = 0, jik = 0; j < N; ++j, ++aj, ++bj) {
                    num_complex xabj = 0.5e0 * (x[aj] + x[bj] + (p[aj] - p[bj]) / (phys::math::im * alpha[j]));
                    num_complex pabj = 0.5e0 * (p[aj] + p[bj] + (phys::math::im * alpha[j]) * (x[aj] - x[bj]));
                    Tab += alpha[j] / (4 * m[j]) + pabj * pabj / (2 * m[j]);
                    for (int i = 0, ik = 0; i < F; ++i) {
                        for (int k = 0; k < F; ++k, ++ik, ++jik) {
                            Vab += std::cobj(cm[i]) * 0.5 * ((xabj - x[bj]) * dVb[jik] + (xabj - x[aj]) * dVa[jik]) *
                                   cn[k];
                        }
                    }
                }
                Hcoeff[ab] = Tab + Vab;
            }
        }
        return 0;
    };

   private:
    num_complex *S, *Sinv;
    num_complex* Snuc* Sele;
    num_complex *S1, *S2;
    num_complex* L;
    num_complex* R;

    num_complex* Hcoeff;
    num_complex* Acoeff;
    num_complex* c;

    void read_param_impl(Param* PM) {
        gamma1 = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));
        gamma2 = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
        xi1    = (1 + Dimension::F * gamma1);
        xi2    = (1 + Dimension::F * gamma2);
    }

    void init_data_impl(DataSet* DS) {
        alpha  = DS->reg<num_real>("integrator.alpha", Dimension::N);
        Acoeff = DS->reg<num_complex>("integrator.Acoeff", Dimension::P);
        Hcoeff = DS->reg<num_complex>("integrator.Hcoeff", Dimension::PP);
        rhored = DS->reg<num_complex>("integrator.rhored", Dimension::FF);
    }

    void init_calc_impl(int stat = -1) {
        for (int iP = 0; iP < Dimension::P; ++iP) {
            //
        }
        exec_kernel(stat);
    }

    int exec_kernel_impl(int stat) {
        calc_Snuc(Snuc, x, p, alpha, Dimension::P, Dimension::N);
        calc_Sele(Sele, c, Dimension::P, Dimension::F);
        for (int ab = 0; ab < Dimension::PP; ++ab) S[ab] = Snuc[ab] * S_ele[ab];
        calc_invS(invS, S, Dimension::P);
        calc_Hcoeff(Hcoeff, );

        for (int ab = 0; ab < Dimension::PP; ++ab) {
            Heff[ab] = S[ab] * (Hcoeff[ab] / Sele[ab] - phys::math::im * dtlnSnuc[ab]);
        }
        ARRAY_MATMUL(Heff, invS, Heff);

        for (int a = 0, ab = 0; a < Dimension::P; ++a) {
            for (int b = 0; b < Dimension::P; ++b, ++ab) { Acoeff[a] += -phys::math::im * Heff[ab] * Acoeff[b] * dt; }
        }
        return 0;
    }
};

};  // namespace PROJECT_NS
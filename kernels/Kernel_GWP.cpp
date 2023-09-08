#include "Kernel_GWP.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec.h"
#include "Kernel_Random.h"
#include "Kernel_Representation.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

namespace PROJECT_NS {

int Kernel_GWP::calc_ekin(num_real* ekin,  // [P]
                          num_real* p,     // [P,N]
                          num_real* m,     // [P,N]
                          int P, int N) {
    Eigen::Map<EigMX<num_real>> Map_ekin(ekin, P, 1);
    Eigen::Map<EigMX<num_real>> Map_p(p, P, N);
    Eigen::Map<EigMX<num_real>> Map_m(m, P, N);
    Map_ekin = (Map_p.array() * Map_p.array() / (2.0e0 * Map_m.array())).matrix().rowwise().sum();
    return 0;
}

/**
 * the expression is exp(-0.25*a*(x1-x2)^2 -0.25*(p1-p2)/a + 0.5i*(p1+p2)(x1-x2) - i(g1-g2))
 */
int Kernel_GWP::calc_Snuc(num_complex* Snuc,  // [P,P]
                          num_real* x1,       // [P,N]
                          num_real* p1,       // [P,N]
                          num_real* m1,       // [P,N]
                          num_real* g1,       // [P]
                          num_real* x2,       // [P,N]
                          num_real* p2,       // [P,N]
                          num_real* m2,       // [P,N]
                          num_real* g2,       // [P]
                          num_real* alpha,    // [N]
                          int P, int N) {
    Eigen::Map<EigMX<num_complex>> Map_Snuc(Snuc, P, P);
    Eigen::Map<EigMX<num_real>> Map_x1(x1, P, N);
    Eigen::Map<EigMX<num_real>> Map_p1(p1, P, N);
    Eigen::Map<EigMX<num_real>> Map_m1(m1, P, N);
    Eigen::Map<EigMX<num_real>> Map_g1(g1, P, 1);

    Eigen::Map<EigMX<num_real>> Map_x2(x2, P, N);
    Eigen::Map<EigMX<num_real>> Map_p2(p2, P, N);
    Eigen::Map<EigMX<num_real>> Map_m2(m2, P, N);
    Eigen::Map<EigMX<num_real>> Map_g2(g2, P, 1);

    Eigen::Map<EigMX<num_real>> Map_a(alpha, N, 1);

    auto O_NP = EigMX<num_real>::Ones(N, P);
    auto O_P1 = EigMX<num_real>::Ones(P, 1);

    auto x1_mul = Map_x1 * Map_a.asDiagonal();
    auto x2_mul = Map_x2 * Map_a.asDiagonal();

    auto p1_div = Map_p1 * (1.0e0 / Map_a.array()).matrix().asDiagonal();
    auto p2_div = Map_p2 * (1.0e0 / Map_a.array()).matrix().asDiagonal();

    auto ax1x1_I = (Map_x1.array() * x1_mul.array()).matrix() * O_NP;
    auto ax1_x2  = x1_mul * Map_x2.transpose();
    auto I_ax2x2 = ((Map_x2.array() * x2_mul.array()).matrix() * O_NP).transpose();
    auto term1   = ax1x1_I + I_ax2x2 - 2 * ax1_x2;

    auto bp1p1_I = (Map_p1.array() * p1_div.array()).matrix() * O_NP;
    auto bp1_p2  = p1_div * Map_p2.transpose();
    auto I_bp2p2 = ((Map_p2.array() * p2_div.array()).matrix() * O_NP).transpose();
    auto term2   = bp1p1_I + I_bp2p2 - 2 * bp1_p2;

    auto x1p1_I = (Map_x1.array() * Map_p1.array()).matrix() * O_NP;
    auto I_x2p2 = ((Map_x2.array() * Map_p2.array()).matrix() * O_NP).transpose();
    auto x1_p2  = Map_x1 * Map_p2.transpose();
    auto p1_x2  = Map_p1 * Map_x2.transpose();
    auto term3  = x1p1_I + x1_p2 - p1_x2 - I_x2p2;

    auto g1_I  = Map_g1 * O_P1.transpose();
    auto I_g2  = O_P1 * Map_g2.transpose();
    auto term4 = g1_I - I_g2;

    Map_Snuc = (-0.25 * (term1 + term2) + 0.5 * phys::math::im * term3 - phys::math::im * term4).array().exp().matrix();

    // // ARRAY_SHOW(Snuc, P, P);
    // for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
    //     for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
    //         num_complex term = 0.0e0;
    //         for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
    //             term += -0.25 * alpha[j] * (x1[aj] - x2[bj]) * (x1[aj] - x2[bj])         //
    //                     - 0.25 / alpha[j] * (p1[aj] - p2[bj]) * (p1[aj] - p2[bj])        //
    //                     + 0.5 * phys::math::im * (p1[aj] + p2[bj]) * (x1[aj] - x2[bj]);  //
    //         }
    //         Snuc[ab] = std::exp(term - phys::math::im * (g1[a] - g2[b]));
    //     }
    // }
    // // ARRAY_SHOW(Snuc, P, P);
    return 0;
}

int Kernel_GWP::calc_Sele(num_complex* Sele,  // [P,P]
                          num_complex* c1,    // [P,F]
                          num_complex* c2,    // [P,F]
                          num_real xi,        //
                          num_real gamma,     //
                          int P, int F) {
    Eigen::Map<EigMX<num_complex>> Map_S(Sele, P, P);
    Eigen::Map<EigMX<num_complex>> Map_c1(c1, P, F);
    Eigen::Map<EigMX<num_complex>> Map_c2(c2, P, F);
    Map_S = xi * Map_c1.conjugate() * Map_c2.transpose() - F * gamma * EigMXc::Identity(P, P);  // or Ones(P, P) @BAD!!!
    return 0;
}

int Kernel_GWP::calc_dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                              num_real* x,            // [P,N]
                              num_real* p,            // [P,N]
                              num_real* m,            // [P,N]
                              num_real* f,            // [P,N]
                              num_real* alpha,        // [N]
                              num_real* ekin,         // [P]
                              int P, int N) {
    Eigen::Map<EigMX<num_complex>> Map_dtlnSnuc(dtlnSnuc, P, P);
    Eigen::Map<EigMX<num_real>> Map_x(x, P, N);
    Eigen::Map<EigMX<num_real>> Map_p(p, P, N);
    Eigen::Map<EigMX<num_real>> Map_m(m, P, N);
    Eigen::Map<EigMX<num_real>> Map_f(f, P, N);
    Eigen::Map<EigMX<num_real>> Map_ekin(ekin, P, 1);
    Eigen::Map<EigMX<num_real>> Map_a(alpha, N, 1);

    auto O_NP = EigMX<num_real>::Ones(N, P);
    auto O_P1 = EigMX<num_real>::Ones(P, 1);

    auto x_mul = Map_x * Map_a.asDiagonal();
    auto p_div = Map_p * (1.0e0 / Map_a.array()).matrix().asDiagonal();

    auto ve     = (Map_p.array() / Map_m.array()).matrix();
    auto ax_ve  = x_mul * ve.transpose();
    auto I_axve = O_NP.transpose() * (x_mul.array() * ve.array()).matrix().transpose();
    auto p_ve   = Map_p * ve.transpose();
    auto I_pve  = O_NP.transpose() * (Map_p.array() * ve.array()).matrix().transpose();

    auto bp_f  = p_div * Map_f.transpose();
    auto I_bpf = O_NP.transpose() * (p_div.array() * Map_f.array()).matrix().transpose();
    auto x_f   = Map_x * Map_f.transpose();
    auto I_xf  = O_NP.transpose() * (Map_x.array() * Map_f.array()).matrix().transpose();

    auto I_ekin = O_P1 * Map_ekin.transpose();

    auto term1 = 0.5e0 * (ax_ve - I_axve) - 0.5e0 * phys::math::im * (p_ve + I_pve);
    auto term2 = -(0.5e0 * (bp_f - I_bpf) + 0.5e0 * phys::math::im * (x_f - I_xf));
    auto term3 = phys::math::im * I_ekin;

    Map_dtlnSnuc = term1 + term2 + term3;

    // // ARRAY_SHOW(dtlnSnuc, P, P);
    // for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
    //     for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
    //         num_complex term = 0.0e0;
    //         for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
    //             term +=
    //                 p[bj] / m[bj] * 0.5 * alpha[j] * (x[aj] - x[bj] + (p[aj] + p[bj]) / (phys::math::im * alpha[j]));
    //             term += (-f[bj]) * 0.5 / alpha[j] * ((p[aj] - p[bj]) + phys::math::im * alpha[j] * (x[aj] - x[bj]));
    //             term += phys::math::im * p[bj] * p[bj] / (2 * m[bj]);
    //         }
    //         dtlnSnuc[ab] = term;
    //     }
    // }
    // // ARRAY_SHOW(dtlnSnuc, P, P);
    return 0;
}

int Kernel_GWP::calc_dtSele(num_complex* dtSele,  // [P,P]
                            num_complex* Sele,    // [P,P]
                            num_complex* c,       // [P,F]
                            num_complex* H,       // [P,F,F]
                            num_real* vpes,       // [P]
                            int P, int F) {
    int FF = F * F;

    // Eigen::Map<EigMX<num_complex>> Map_dtSele(dtSele, P, P);
    // Eigen::Map<EigMX<num_complex>> Map_S(Sele, P, P);
    // Eigen::Map<EigMX<num_complex>> Map_c(c, P, F);
    // Eigen::Map<EigMX<num_complex>> Map_H_P_FF(H, P, F * F);
    // Eigen::Map<EigMX<num_complex>> Map_H_PF_F(H, P * F, F);
    // Eigen::Map<EigMX<num_complex>> Map_H_F_PF(H, F, P * F);
    // Eigen::Map<EigMX<num_real>> Map_vpes(vpes, P, 1);

    // Map_S = Map_c.conjugate() * Map_c.transpose();
    for (int a = 0, ab = 0, aF = 0; a < P; ++a, aF += F) {
        num_complex* ca = c + aF;
        for (int b = 0, bF = 0, bFF = 0; b < P; ++b, ++ab, bF += F, bFF += FF) {
            num_complex* cb = c + bF;
            num_complex* Hb = H + bFF;

            num_complex term1 = ARRAY_INNER_VMV_TRANS1(ca, H, cb, F, F);
            num_complex term2 = ARRAY_INNER_TRANS1(ca, cb, F);
            dtSele[ab]        = -phys::math::im * (term1 + vpes[b] * term2);
        }
    }
    return 0;
}

int Kernel_GWP::calc_invS(num_complex* invS, num_complex* S, int P) {
    Eigen::Map<EigMX<num_complex>> Map_invS(invS, Dimension::P, Dimension::P);
    Eigen::Map<EigMX<num_complex>> Map_S(S, Dimension::P, Dimension::P);

    Map_invS = Map_S.lu().inverse();

    return 0;
}

double Kernel_GWP::calc_density(num_complex* rhored,         // [F,F]
                                num_complex* Acoeff,         // [P]
                                num_complex* Snuc,           // [P,P]
                                num_complex* c,              // [P,F]
                                num_real xi,                 // for kernel
                                num_real gamma,              // for kernel
                                int P_used, int P, int F) {  //@bug
    Eigen::Map<EigMX<num_complex>> Map_rhored(rhored, Dimension::F, Dimension::F);
    Eigen::Map<EigMX<num_complex>> Map_Acoeff(Acoeff, Dimension::P, 1);
    Eigen::Map<EigMX<num_complex>> Map_Snuc(Snuc, Dimension::P, Dimension::P);
    Eigen::Map<EigMX<num_complex>> Map_c(c, Dimension::P, Dimension::F);
    Map_rhored =
        (xi * Map_c.adjoint() * Map_Acoeff.adjoint().asDiagonal() * Map_Snuc * Map_Acoeff.asDiagonal() * Map_c -
         gamma * (Map_Acoeff.adjoint() * Map_Snuc * Map_Acoeff).sum() * EigMXc::Identity(F, F));
    return std::abs(Map_rhored.trace());
}

int Kernel_GWP::calc_Hbasis(num_complex* Hbasis,  // [P,P]
                            num_real* vpes,       // [P]
                            num_real* grad,       // [P,N]
                            num_real* V,          // [P,F,F]
                            num_real* dV,         // [P,N,F,F]
                            num_real* x,          // [P,N]
                            num_real* p,          // [P,N]
                            num_real* m,          // [P,N]
                            num_real* alpha,      // [N]
                            num_complex* Sele,    // [P,P]
                            num_complex* c,       // [P,F]
                            int P, int N, int F   // Dimensions
) {
    int FF  = F * F;
    int NFF = N * FF;

    for (int a = 0, aN = 0, aF = 0, ab = 0; a < P; ++a, aN += N, aF += F) {
        num_real* Va    = V + a * FF;
        num_complex* ca = c + aF;
        for (int b = 0, bN = 0, bF = 0; b < P; ++b, ++ab, bN += N, bF += F) {
            num_real* Vb    = V + b * FF;
            num_complex* cb = c + bF;

            num_complex Tab = 0.0e0;
            num_complex Vab = 0.0e0;

            Vab += 0.5e0 * (vpes[a] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca, Va, cb, F, F));
            Vab += 0.5e0 * (vpes[b] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca, Vb, cb, F, F));
            for (int j = 0, jFF = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj, jFF += FF) {
                num_real* dVaj   = dV + aj * FF;
                num_real* dVbj   = dV + bj * FF;
                num_complex xabj = 0.5e0 * (x[aj] + x[bj] + (p[aj] - p[bj]) / (phys::math::im * alpha[j]));
                num_complex pabj = 0.5e0 * (p[aj] + p[bj] + (phys::math::im * alpha[j]) * (x[aj] - x[bj]));
                Tab += alpha[j] / (4 * m[bj]) + pabj * pabj / (2 * m[bj]);

                Vab += 0.5e0 * (xabj - x[aj]) * (grad[aj] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca, dVaj, cb, F, F));
                Vab += 0.5e0 * (xabj - x[bj]) * (grad[bj] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca, dVbj, cb, F, F));
            }
            Hbasis[ab] = Tab + Vab;
        }
    }
    return 0;
};

int Kernel_GWP::calc_Hbasis_adia(num_complex* Hbasis,  // [P,P]
                                 num_real* E,          // [P,F]
                                 num_real* dE,         // [P,N,F,F]
                                 num_real* x,          // [P,N]
                                 num_real* p,          // [P,N]
                                 num_real* m,          // [P,N]
                                 num_real* alpha,      // [N]
                                 num_complex* c,       // [P,F]
                                 int P, int N, int F   // Dimensions
) {
    int FF  = F * F;
    int NFF = N * FF;

    for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
            num_real* Ea    = E + a * F;
            num_real* Eb    = E + b * F;
            num_real* dEa   = dE + a * NFF;
            num_real* dEb   = dE + b * NFF;
            num_complex* ca = c + a * F;
            num_complex* cb = c + b * F;

            num_complex Tab = 0.0e0;
            num_complex Eab = 0.0e0;

            for (int j = 0, jik = 0, jFF = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj, jFF += Dimension::FF) {
                num_real* dEaj   = dEa + jFF;
                num_real* dEbj   = dEb + jFF;
                num_complex xabj = 0.5e0 * (x[aj] + x[bj] + (p[aj] - p[bj]) / (phys::math::im * alpha[j]));
                num_complex pabj = 0.5e0 * (p[aj] + p[bj] + (phys::math::im * alpha[j]) * (x[aj] - x[bj]));
                Tab += alpha[j] / (4 * m[j]) + pabj * pabj / (2 * m[j]);

                for (int i = 0, ik = 0; i < F; ++i) {
                    for (int k = 0; k < F; ++k, ++ik) {
                        Eab += (i == k) ? 0.5e0 * std::conj(ca[i]) * cb[k] *
                                              (Ea[i] + Eb[i] + (xabj - x[aj]) * dEaj[ik] + (xabj - x[bj]) * dEbj[ik])
                                        : 0.5e0 * phys::math::im * std::conj(ca[i]) * cb[k] *
                                              (p[aj] / m[aj] * dEaj[ik] / (Ea[k] - Ea[i]) +
                                               p[bj] / m[bj] * dEbj[ik] / (Eb[k] - Eb[i]));
                    }
                }
            }
            Hbasis[ab] = Tab + Eab;
        }
    }
    return 0;
}

void Kernel_GWP::read_param_impl(Param* PM) {
    dt            = PM->get<double>("dt", LOC(), phys::time_d);     //
    alpha0        = PM->get<double>("alpha0", LOC(), 1.0f);         //
    width_scaling = PM->get<double>("width_scaling", LOC(), 1.0f);  //
    break_thres   = PM->get<double>("break_thres", LOC(), 1.0f);    //
    P_used0       = PM->get<int>("P_initial", LOC(), 1);            //
    max_clone     = PM->get<int>("max_clone", LOC(), 0);            //
    gamma         = PM->get<double>("gamma", LOC(), 0.0f);          //
    xi            = 1 + Dimension::F * gamma;

    impl_type          = PM->get<int>("impl_type", LOC(), 0);           //
    samp_type          = PM->get<int>("samp_type", LOC(), 0);           //
    aset_type          = PM->get<int>("aset_type", LOC(), 0);           //
    time_displace_step = PM->get<int>("time_displace_step", LOC(), 1);  //
}

void Kernel_GWP::init_data_impl(DataSet* DS) {
    alpha    = DS->reg<num_real>("integrator.alpha", Dimension::N);
    Xcoeff   = DS->reg<num_complex>("integrator.Xcoeff", Dimension::P);
    Acoeff   = DS->reg<num_complex>("integrator.Acoeff", Dimension::P);
    dtAcoeff = DS->reg<num_complex>("integrator.dtAcoeff", Dimension::P);
    Hcoeff   = DS->reg<num_complex>("integrator.Hcoeff", Dimension::PP);
    Hbasis   = DS->reg<num_complex>("integrator.Hbasis", Dimension::PP);
    UXdt     = DS->reg<num_complex>("integrator.UXdt", Dimension::PP);
    rhored   = DS->reg<num_complex>("integrator.rhored", Dimension::FF);
    rhored2  = DS->reg<num_complex>("integrator.rhored2", Dimension::FF);
    rhored3  = DS->reg<num_complex>("integrator.rhored3", Dimension::FF);

    Snuc     = DS->reg<num_complex>("integrator.Snuc", Dimension::PP);
    Sele     = DS->reg<num_complex>("integrator.Sele", Dimension::PP);
    S        = DS->reg<num_complex>("integrator.S", Dimension::PP);
    invS     = DS->reg<num_complex>("integrator.invS", Dimension::PP);
    dtlnSnuc = DS->reg<num_complex>("integrator.dtlnSnuc", Dimension::PP);
    dtSele   = DS->reg<num_complex>("integrator.dtSele", Dimension::PP);
    L        = DS->reg<num_real>("integrator.GWP.L", Dimension::P);
    L1       = DS->reg<num_real>("integrator.GWP.L1", Dimension::P);
    L2       = DS->reg<num_real>("integrator.GWP.L2", Dimension::P);
    R        = DS->reg<num_complex>("integrator.GWP.R", Dimension::PP);
    R1       = DS->reg<num_complex>("integrator.GWP.R1", Dimension::PP);
    R2       = DS->reg<num_complex>("integrator.GWP.R2", Dimension::PP);
    S1       = DS->reg<num_complex>("integrator.GWP.S1", Dimension::PP);
    S1h      = DS->reg<num_complex>("integrator.GWP.S1h", Dimension::PP);
    invS1h   = DS->reg<num_complex>("integrator.GWP.invS1h", Dimension::PP);
    S2       = DS->reg<num_complex>("integrator.GWP.S2", Dimension::PP);
    S2h      = DS->reg<num_complex>("integrator.GWP.S2h", Dimension::PP);
    invS2h   = DS->reg<num_complex>("integrator.GWP.invS2h", Dimension::PP);
    Sx       = DS->reg<num_complex>("integrator.GWP.Sx", Dimension::PP);

    ekin          = DS->reg<num_real>("integrator.ekin", Dimension::P);
    g             = DS->reg<num_real>("integrator.g", Dimension::P);
    clone_account = DS->reg<int>("integrator.clone_account", Dimension::P);
    // pf_cross      = DS->reg<bool>("integrator.pf_cross", Dimension::PF);

    //
    Udt  = DS->reg<num_complex>("integrator.Udt", Dimension::PFF);
    H    = DS->reg<num_complex>("model.rep.H", Dimension::PFF);
    vpes = DS->reg<num_real>("model.vpes", Dimension::P);
    grad = DS->reg<num_real>("model.grad", Dimension::PN);
    V    = DS->reg<num_real>("model.V", Dimension::PFF);
    dV   = DS->reg<num_real>("model.dV", Dimension::PNFF);
    E    = DS->reg<num_real>("model.rep.E", Dimension::PF);
    T    = DS->reg<num_real>("model.rep.T", Dimension::PFF);
    dE   = DS->reg<num_real>("model.rep.dE", Dimension::PNFF);
    x    = DS->reg<num_real>("integrator.x", Dimension::PN);
    p    = DS->reg<num_real>("integrator.p", Dimension::PN);
    m    = DS->reg<num_real>("integrator.m", Dimension::PN);
    f    = DS->reg<num_real>("integrator.f", Dimension::PN);
    c    = DS->reg<num_complex>("integrator.c", Dimension::PF);

    fun_diag_F = DS->reg<num_complex>("integrator.tmp.fun_diag_F", Dimension::F);
    fun_diag_P = DS->reg<num_complex>("integrator.tmp.fun_diag_P", Dimension::P);
    MatR_PP    = DS->reg<num_real>("integrator.tmp.MatR_PP", Dimension::PP);
    MatC_PP    = DS->reg<num_complex>("integrator.tmp.MatC_PP", Dimension::PP);
    I_PP       = DS->reg<num_complex>("integrator.tmp.I_PP", Dimension::PP);
    Ubranch    = DS->reg<num_complex>("integrator.Ubranch", Dimension::FF);

    x_last    = DS->reg<num_real>("integrator.x_last", Dimension::PN);
    p_last    = DS->reg<num_real>("integrator.p_last", Dimension::PN);
    grad_last = DS->reg<num_real>("integrator.grad_last", Dimension::PN);
    dV_last   = DS->reg<num_real>("integrator.dV_last", Dimension::PNFF);
    g_last    = DS->reg<num_real>("integrator.g_last", Dimension::P);
    c_last    = DS->reg<num_complex>("integrator.c_last", Dimension::PF);

    P_used_ptr = DS->reg<num_real>("integrator.P_used");
    norm_ptr   = DS->reg<num_real>("integrator.norm");
    gammat_ptr = DS->reg<num_real>("integrator.gammat");
    veF        = DS->reg<num_real>("integrator.veF", Dimension::FF);
}

void Kernel_GWP::init_calc_impl(int stat) {
    // switch (samp_type) {
    //     case 0: {
    //         break;
    //     }
    //     case 1:
    //     case 2: {
    //         break;
    //     }
    //     case 3: {
    //         break;
    //     }
    // }
    if (samp_type < 3) {  // overlap or neighbourhood re-sampling
        for (int iP = 0; iP < Dimension::P; ++iP) {
            num_complex* w = Kernel_Elec::w + iP;
            num_complex* c = Kernel_Elec::c + iP * Dimension::F;
            num_complex* U = Kernel_Elec::U + iP * Dimension::FF;
            int* occ_nuc   = Kernel_Elec::occ_nuc + iP;

            num_real* x = this->x + iP * Dimension::N;
            num_real* p = this->p + iP * Dimension::N;

            /////////////////////////////////////////////////////////////////
            if (samp_type == 1)
                for (int j = 0; j < Dimension::N; ++j) {
                    x[j] = this->x[j];
                    p[j] = this->p[j];
                }

            if (samp_type == 2 && iP > 0) {
                for (int j = 0; j < Dimension::N; ++j) {
                    double randu;
                    Kernel_Random::rand_gaussian(&randu);
                    x[j] = this->x[j] + width_scaling * randu / sqrt(Dimension::N * alpha0);
                    Kernel_Random::rand_gaussian(&randu);
                    p[j] = this->p[j] + width_scaling * randu * sqrt(alpha0 / Dimension::N);
                }
            }

            *w       = 1.0e0;              ///< initial measure
            *occ_nuc = Kernel_Elec::occ0;  ///< initial occupation
            /// < initial c (not used)
            for (int i = 0; i < Dimension::F; ++i) {
                double randu;
                Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
                c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
            }
            ARRAY_EYE(U, Dimension::F);  ///< initial propagator
        }
    }

    /// time displaced re-sampling
    if (samp_type == 3) {
        *Kernel_Elec::w       = 1.0e0;
        *Kernel_Elec::occ_nuc = Kernel_Elec::occ0;
        for (int i = 0; i < Dimension::F; ++i) {
            double randu;
            Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
            Kernel_Elec::c[i] =
                (i == *Kernel_Elec::occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
        }
        ARRAY_EYE(Kernel_Elec::U, Dimension::F);
        for (int iP = 1; iP < Dimension::P; ++iP) {
            num_real* x_now          = x + iP * Dimension::N;
            num_real* p_now          = p + iP * Dimension::N;
            num_real* f_now          = f + iP * Dimension::N;
            num_complex* U_now       = Kernel_Elec::U + iP * Dimension::FF;
            num_complex* c_now       = Kernel_Elec::c + iP * Dimension::F;
            num_complex* rho_nuc_now = Kernel_Elec::rho_nuc + iP * Dimension::FF;

            num_real* x_prev    = x + std::max({iP - 2, 0}) * Dimension::N;
            num_real* p_prev    = p + std::max({iP - 2, 0}) * Dimension::N;
            num_real* f_prev    = f + std::max({iP - 2, 0}) * Dimension::N;
            num_complex* U_prev = Kernel_Elec::U + std::max({iP - 2, 0}) * Dimension::FF;

            num_real* E_now      = E + iP * Dimension::F;
            num_real* T_now      = T + iP * Dimension::FF;
            num_complex* Udt_now = Udt + iP * Dimension::FF;

            num_real signdt = (iP % 2 == 0) ? dt : -dt;

            for (int j = 0; j < Dimension::N; ++j) x_now[j] = x_prev[j], p_now[j] = p_prev[j], f_now[j] = f_prev[j];

            for (int istep_displace = 0; istep_displace < time_displace_step; ++istep_displace) {
                for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
                for (int j = 0; j < Dimension::N; ++j) x_now[j] += p_now[j] / m[j] * signdt;
                _kmodel->exec_kernel();  // only iP needed
                _krepr->exec_kernel();   // only iP needed
                switch (Kernel_Representation::ele_repr_type) {
                    case RepresentationPolicy::Diabatic: {
                        for (int i = 0; i < Dimension::F; ++i) fun_diag_F[i] = exp(-phys::math::im * E_now[i] * signdt);
                        ARRAY_MATMUL3_TRANS2(Udt_now, T_now, fun_diag_F, T_now, Dimension::F, Dimension::F, 0,
                                             Dimension::F);
                        break;
                    }
                }
                ARRAY_MATMUL(U_now, Udt_now, U_prev, Dimension::F, Dimension::F, Dimension::F);
                ARRAY_MATMUL(c_now, U_now, c, Dimension::F, Dimension::F, 1);
                Kernel_Elec::ker_from_c(rho_nuc_now, c_now, xi, gamma, Dimension::F);
                _kforce->exec_kernel();
                for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
            }
            // ARRAY_SHOW(x_now, 1, Dimension::N);
            // ARRAY_SHOW(p_now, 1, Dimension::N);
        }

        for (int iP = 0; iP < Dimension::P; ++iP) {
            num_complex* w = Kernel_Elec::w + iP;
            num_complex* c = Kernel_Elec::c + iP * Dimension::F;
            num_complex* U = Kernel_Elec::U + iP * Dimension::FF;
            int* occ_nuc   = Kernel_Elec::occ_nuc + iP;

            /////////////////////////////////////////////////////////////////

            *w       = 1.0e0;              ///< initial measure
            *occ_nuc = Kernel_Elec::occ0;  ///< initial occupation
            /// < use time-displaced c ?
            // for (int i = 0; i < Dimension::F; ++i) {
            //     double randu;
            //     Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
            //     c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
            // }
            ARRAY_EYE(U, Dimension::F);  ///< initial propagator reset to identity
        }
    }

    for (int i = 0; i < Dimension::N; ++i) alpha[i] = alpha0;
    for (int a = 0; a < Dimension::P; ++a) g[a] = 0.0e0;
    if (aset_type == 0) {
        for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = (a == 0) ? 1.0e0 : 0.0e0;
    } else if (aset_type == 1) {
        for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = 1.0e0;  // @only for test
    }

    // normalization of A
    P_used        = P_used0;
    P_used_ptr[0] = P_used;
    // std::cout << "P_used_INIT = " << P_used << "\n";
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) { S[ab] = Snuc[ab] * Sele[ab]; }
    num_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff, S, Acoeff, 1, Dimension::P, Dimension::P, 1);
    for (int a = 0; a < Dimension::P; ++a) Acoeff[a] /= sqrt(abs(scale));
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff, Snuc, Acoeff, 1, Dimension::P, Dimension::P, 1);
    xi = 1.0e0 + Dimension::F * gamma * abs(scale);

    norm_ptr[0] = 1.0e0;

    // ARRAY_SHOW(Acoeff, 1, Dimension::P);

    Kernel_Elec::c_init = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::T_init = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);

    for (int a = 0; a < Dimension::P; ++a) g_last[a] = g[a];
    for (int aj = 0; aj < Dimension::PN; ++aj) x_last[aj] = x[aj];
    for (int aj = 0; aj < Dimension::PN; ++aj) p_last[aj] = p[aj];
    for (int ai = 0; ai < Dimension::PF; ++ai) c_last[ai] = c[ai];

    double dt_backup = dt;
    dt               = 0.0e0;
    exec_kernel();  // don't evolve!!!
    dt = dt_backup;
}

int Kernel_GWP::impl_0(int stat) {
    // calculate overlap integral & time derivative of overlap integral

    // ARRAY_SHOW(p, Dimension::P, Dimension::N);
    // ARRAY_SHOW(x, Dimension::P, Dimension::N);
    // ARRAY_SHOW(c, Dimension::P, Dimension::F);
    // ARRAY_SHOW(g, 1, Dimension::P);

    calc_ekin(ekin, p, m, Dimension::P, Dimension::N);
    calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
    calc_dtlnSnuc(dtlnSnuc, x, p, m, f, alpha, ekin, Dimension::P, Dimension::N);
    calc_dtSele(dtSele, Sele, c, H, vpes, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) { S[ab] = Snuc[ab] * Sele[ab]; }
    // calculate inverse of overlap integral
    EigenSolve(L1, R1, S, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / L1[a];
    ARRAY_MATMUL3_TRANS2(invS, R1, fun_diag_P, R1, Dimension::P, Dimension::P, 0, Dimension::P);

    // ARRAY_SHOW(Snuc, Dimension::P, Dimension::P);
    // ARRAY_SHOW(Sele, Dimension::P, Dimension::P);
    // ARRAY_SHOW(S, Dimension::P, Dimension::P);
    // ARRAY_SHOW(invS, Dimension::P, Dimension::P);

    // calculate Hamiltonian between two configurations basis
    calc_Hbasis(Hbasis, vpes, grad, V, dV, x, p, m, alpha, Sele, c, Dimension::P, Dimension::N, Dimension::F);
    // calculate effective Hamiltonian
    for (int ab = 0; ab < Dimension::PP; ++ab) {
        Hcoeff[ab] =
            (Snuc[ab] * Hbasis[ab] - phys::math::im * S[ab] * dtlnSnuc[ab] - phys::math::im * Snuc[ab] * dtSele[ab]);
    }
    ARRAY_MATMUL(Hcoeff, invS, Hcoeff, Dimension::P, Dimension::P, Dimension::P);

    // ARRAY_SHOW(Hcoeff, Dimension::P, Dimension::P);

    Eigen::Map<EigMXc> Map_Hcoeff(Hcoeff, Dimension::P, Dimension::P);
    Eigen::Map<EigMXc> Map_UXdt(UXdt, Dimension::P, Dimension::P);
    auto eigr = Eigen::ComplexEigenSolver<EigMXc>(Map_Hcoeff);
    auto Er   = eigr.eigenvalues();
    auto Vr   = eigr.eigenvectors();
    auto eigl = Eigen::ComplexEigenSolver<EigMXc>(Map_Hcoeff.adjoint());
    auto El   = eigl.eigenvalues();
    auto Vl   = eigl.eigenvectors();
    auto Slr  = (Vl.adjoint() * Vr).diagonal();
    Map_UXdt  = Vr * ((-phys::math::im * Er.array() * dt).exp() / Slr.array()).matrix().asDiagonal() * Vl.adjoint();

    // ARRAY_SHOW(UXdt, Dimension::P, Dimension::P);

    // update Acoeff
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    ARRAY_MATMUL(Acoeff, UXdt, Acoeff, Dimension::P, Dimension::P, 1);
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    num_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff, S, Acoeff, 1, Dimension::P, Dimension::P, 1);
    norm_ptr[0] *= abs(scale);
    for (int a = 0; a < Dimension::P; ++a) Acoeff[a] /= sqrt(abs(scale));

    for (int a = 0; a < Dimension::P; ++a) g[a] += ekin[a] * dt;

    cloning();
    death();
    return 0;
}

int Kernel_GWP::impl_1(int stat) {
    // ARRAY_SHOW(x_last, Dimension::P, Dimension::N);
    // ARRAY_SHOW(x, Dimension::P, Dimension::N);

    // calculate overlap integral
    calc_ekin(ekin, p, m, Dimension::P, Dimension::N);
    calc_Snuc(Snuc, x_last, p_last, m, g_last, x_last, p_last, m, g_last, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c_last, c_last, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) S1[ab] = Snuc[ab] * Sele[ab];

    // ARRAY_SHOW(S1, Dimension::P, Dimension::P);

    calc_Snuc(Snuc, x, p, m, g, x_last, p_last, m, g_last, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c_last, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) Sx[ab] = Snuc[ab] * Sele[ab];

    // ARRAY_SHOW(Sx, Dimension::P, Dimension::P);

    calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) S2[ab] = Snuc[ab] * Sele[ab];
    for (int ab = 0; ab < Dimension::PP; ++ab) S[ab] = S2[ab];

    // ARRAY_SHOW(S2, Dimension::P, Dimension::P);

    EigenSolve(L1, R1, S1, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = sqrt(L1[a]);
    ARRAY_MATMUL3_TRANS2(S1h, R1, fun_diag_P, R1, Dimension::P, Dimension::P, 0, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / sqrt(L1[a]);
    ARRAY_MATMUL3_TRANS2(invS1h, R1, fun_diag_P, R1, Dimension::P, Dimension::P, 0, Dimension::P);

    EigenSolve(L2, R2, S2, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = sqrt(L2[a]);
    ARRAY_MATMUL3_TRANS2(S2h, R2, fun_diag_P, R2, Dimension::P, Dimension::P, 0, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / sqrt(L2[a]);
    ARRAY_MATMUL3_TRANS2(invS2h, R2, fun_diag_P, R2, Dimension::P, Dimension::P, 0, Dimension::P);

    ARRAY_MATMUL(MatC_PP, invS2h, invS2h, Dimension::P, Dimension::P, Dimension::P);
    ARRAY_MATMUL3_TRANS1(MatC_PP, Sx, MatC_PP, Sx, Dimension::P, Dimension::P, Dimension::P, Dimension::P);
    double error = 0.0e0;
    for (int ab = 0; ab < Dimension::PP; ++ab) MatC_PP[ab] -= S1[ab], error += abs(MatC_PP[ab]);
    // ARRAY_SHOW(MatC_PP, Dimension::P, Dimension::P);
    std::cout << "test_dt2: " << error << "\n";

    // calculate Hamiltonian between two configurations basis
    calc_Hbasis(Hbasis, vpes, grad, V, dV, x, p, m, alpha, Sele, c, Dimension::P, Dimension::N, Dimension::F);
    ARRAY_MATMUL(Hcoeff, invS1h, Hbasis, Dimension::P, Dimension::P, Dimension::P);
    ARRAY_MATMUL(Hcoeff, Hcoeff, invS1h, Dimension::P, Dimension::P, Dimension::P);

    // ARRAY_SHOW(Hcoeff, Dimension::P, Dimension::P);

    EigenSolve(L, R, Hcoeff, Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = exp(-phys::math::im * L[a] * dt);
    ARRAY_MATMUL3_TRANS2(UXdt, R, fun_diag_P, R, Dimension::P, Dimension::P, 0, Dimension::P);

    ARRAY_MATMUL(Xcoeff, S1h, Acoeff, Dimension::P, Dimension::P, 1);

    ARRAY_MATMUL(Xcoeff, UXdt, Xcoeff, Dimension::P, Dimension::P, 1);
    ARRAY_MATMUL(Xcoeff, invS1h, Xcoeff, Dimension::P, Dimension::P, 1);
    ARRAY_MATMUL(Xcoeff, Sx, Xcoeff, Dimension::P, Dimension::P, 1);
    ARRAY_MATMUL(Xcoeff, invS2h, Xcoeff, Dimension::P, Dimension::P, 1);

    ARRAY_MATMUL(Acoeff, invS2h, Xcoeff, Dimension::P, Dimension::P, 1);

    num_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff, S, Acoeff, 1, Dimension::P, Dimension::P, 1);
    norm_ptr[0] *= abs(scale);

    // update phase
    for (int a = 0; a < Dimension::P; ++a) g[a] += ekin[a] * dt;

    cloning();
    death();

    for (int a = 0; a < Dimension::P; ++a) g_last[a] = g[a];
    for (int aj = 0; aj < Dimension::PN; ++aj) x_last[aj] = x[aj];
    for (int aj = 0; aj < Dimension::PN; ++aj) p_last[aj] = p[aj];
    for (int ai = 0; ai < Dimension::PF; ++ai) c_last[ai] = c[ai];
    return 0;
}
int Kernel_GWP::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* c_init  = Kernel_Elec::c_init + iP * Dimension::F;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* K1      = Kernel_Elec::K1 + iP * Dimension::FF;

        num_real* T      = Kernel_Elec::T + iP * Dimension::FF;
        num_real* T_init = Kernel_Elec::T_init + iP * Dimension::FF;

        /////////////////////////////////////////////////////////////////

        for (int i = 0; i < Dimension::F; ++i) c[i] = c_init[i];

        // 1) transform from inp_repr => ele_repr
        Kernel_Representation::transform(c, T_init, Dimension::F,               //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::H);
        // 2) propagte along ele_repr
        if (stat >= 0) ARRAY_MATMUL(c, U, c, Dimension::F, Dimension::F, 1);

        // 3) transform back from ele_repr => inp_repr
        Kernel_Representation::transform(c, T, Dimension::F,                    //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::H);

        Kernel_Elec::ker_from_c(K1, c, xi, gamma, Dimension::F);
        Kernel_Elec::ker_from_c(rho_nuc, c, xi, gamma, Dimension::F);
    }
    switch (impl_type) {
        case 0: {
            impl_0(stat);
            break;
        }
        case 1: {
            impl_1(stat);
            break;
        }
    }
    num_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff, Snuc, Acoeff, 1, Dimension::P, Dimension::P, 1);
    xi = 1.0e0 + Dimension::F * gamma * abs(scale);
    calc_density(rhored, Acoeff, Snuc, c, xi, gamma, P_used, Dimension::P, Dimension::F);

    ARRAY_CLEAR(rhored2, Dimension::FF);
    for (int a = 0, aik = 0; a < P_used; ++a) {
        for (int ik = 0; ik < Dimension::FF; ++ik, ++aik) rhored2[ik] += Kernel_Elec::K1[aik];
    }
    for (int ik = 0; ik < Dimension::FF; ++ik) rhored2[ik] /= double(P_used);
    return 0;
}

int Kernel_GWP::cloning() {
    P_used_ptr[0] = P_used;
    if (P_used >= Dimension::P) return 0;
    int P_increase = P_used;

    // calc_Snuc(Snuc, this->x, this->p, this->m, this->g, this->x, this->p, this->m, this->g, alpha, Dimension::P,
    //           Dimension::N);
    // calc_density(rhored, Acoeff, Snuc, c, xi, gamma, P_used, Dimension::P, Dimension::F);
    // ARRAY_SHOW(Snuc, Dimension::P, Dimension::P);
    // ARRAY_SHOW(c, Dimension::P, Dimension::F);
    // ARRAY_SHOW(Acoeff, Dimension::P, 1);
    // ARRAY_SHOW(rhored, Dimension::F, Dimension::F);

    // calc_density(rhored, Acoeff, Snuc, c, 1, 0, P_used, Dimension::P, Dimension::F);
    // ARRAY_SHOW(rhored, Dimension::F, Dimension::F);

    for (int iP = 0; iP < P_used; ++iP) {
        num_real* g         = this->g + iP;
        num_real* x         = this->x + iP * Dimension::N;
        num_real* p         = this->p + iP * Dimension::N;
        num_real* f         = this->f + iP * Dimension::N;
        num_real* grad      = this->grad + iP * Dimension::N;
        num_real* V         = this->V + iP * Dimension::FF;
        num_real* dV        = this->dV + iP * Dimension::NFF;
        num_complex* c      = this->c + iP * Dimension::F;
        num_complex* c_init = Kernel_Elec::c_init + iP * Dimension::F;
        num_complex* U      = Kernel_Elec::U + iP * Dimension::FF;

        /////////////////////////////////////////////////

        double max_break_val = 0.0e0;
        int break_state_i = -1, break_state_k = -1;
        Eigen::Map<EigMXr> Map_veF(veF, Dimension::FF, 1);
        Eigen::Map<EigMXr> Map_dV(dV, Dimension::N, Dimension::FF);
        Eigen::Map<EigMXr> Map_p(p, Dimension::N, 1);
        Eigen::Map<EigMXr> Map_m(m, Dimension::N, 1);
        Map_veF = Map_dV.transpose() * (Map_p.array() / Map_m.array()).matrix();

        for (int i = 0; i < Dimension::F; ++i) {
            for (int k = i + 1; k < Dimension::F; ++k) {
                double break_val =
                    std::abs(c[i] * c[i] * c[k] * c[k] * (veF[i * Dimension::Fadd1] - veF[k * Dimension::Fadd1]) /
                             (1e-6 + abs(V[i * Dimension::Fadd1] - V[k * Dimension::Fadd1])));
                if (break_val > max_break_val) {
                    max_break_val = break_val;
                    break_state_i = i;
                    break_state_k = k;
                }
            }
        }
        if (max_break_val > break_thres) {
            int break_state = (std::abs(c[break_state_i]) > std::abs(c[break_state_k])) ? break_state_i : break_state_k;
            double norm_b   = abs(c[break_state]);
            double norm_a   = sqrt(1 - norm_b * norm_b);
            num_complex norm_phase_a = norm_a * (norm_a + phys::math::im * norm_b);
            num_complex norm_phase_b = norm_b * (norm_b - phys::math::im * norm_a);

            // std::cout << "norm_a, norm_b, iP, P_used, P_increase: " << norm_a << ", " << norm_b << ", " << iP << ","
            //           << P_used << ", " << P_increase << "\n";

            num_real* g_new         = this->g + P_increase;
            num_real* x_new         = this->x + P_increase * Dimension::N;
            num_real* p_new         = this->p + P_increase * Dimension::N;
            num_real* f_new         = this->f + P_increase * Dimension::N;
            num_real* grad_new      = this->grad + P_increase * Dimension::N;
            num_real* dV_new        = this->dV + P_increase * Dimension::NFF;
            num_complex* c_new      = this->c + P_increase * Dimension::F;
            num_complex* c_init_new = Kernel_Elec::c_init + P_increase * Dimension::F;
            num_complex* U_new      = Kernel_Elec::U + P_increase * Dimension::FF;

            g_new[0] = g[0];
            for (int j = 0; j < Dimension::N; ++j) x_new[j] = x[j];
            for (int j = 0; j < Dimension::N; ++j) p_new[j] = p[j];
            for (int j = 0; j < Dimension::N; ++j) f_new[j] = f[j];
            for (int j = 0; j < Dimension::N; ++j) grad_new[j] = grad[j];
            for (int j = 0; j < Dimension::NFF; ++j) dV_new[j] = dV[j];

            for (int i = 0; i < Dimension::F; ++i) fun_diag_F[i] = c[i];
            for (int i = 0; i < Dimension::F; ++i) c_new[i] = ((i == break_state) ? 0.0e0 : (c[i] / norm_phase_a));
            for (int i = 0; i < Dimension::F; ++i) c_init_new[i] = c_init[i];
            ARRAY_MATMUL_TRANS2(Ubranch, c_new, fun_diag_F, Dimension::F, 1, Dimension::F);
            ARRAY_MATMUL(U_new, Ubranch, U, Dimension::F, Dimension::F, Dimension::F);

            for (int i = 0; i < Dimension::F; ++i) c[i] = ((i == break_state) ? c[i] / norm_phase_b : 0.0e0);
            ARRAY_MATMUL_TRANS2(Ubranch, c, fun_diag_F, Dimension::F, 1, Dimension::F);
            ARRAY_MATMUL(U, Ubranch, U, Dimension::F, Dimension::F, Dimension::F);

            Acoeff[P_increase] = Acoeff[iP] * norm_phase_a;
            Acoeff[iP] *= norm_phase_b;

            P_increase++;

            if (P_increase >= Dimension::P) break;
        }
    }

    // std::cout << "cloning " << P_increase - P_used << " trajs\n";
    // std::cout << "P_used = " << P_used << "\n";

    if (P_increase > P_used) {
        // std::cout << "\n#############################################\n\n";

        P_used = P_increase;
        // calc_Snuc(Snuc, this->x, this->p, this->m, this->g, this->x, this->p, this->m, this->g, alpha, Dimension::P,
        //           Dimension::N);
        // calc_density(rhored, Acoeff, Snuc, c, xi, gamma, P_used, Dimension::P, Dimension::F);
        // ARRAY_SHOW(Snuc, Dimension::P, Dimension::P);
        // ARRAY_SHOW(c, Dimension::P, Dimension::F);
        // ARRAY_SHOW(Acoeff, Dimension::P, 1);
        // ARRAY_SHOW(rhored, Dimension::F, Dimension::F);
        // calc_density(rhored, Acoeff, Snuc, c, 1, 0, P_used, Dimension::P, Dimension::F);
        // ARRAY_SHOW(rhored, Dimension::F, Dimension::F);
    }


    // std::cout << "after P_used = " << P_used << "\n";
    return 0;

    // for (int iP = 0; iP < Dimension::P; ++iP) {
    //     num_real* p    = this->p + iP * Dimension::N;
    //     num_real* grad = this->grad + iP * Dimension::N;
    //     num_real* dV   = this->dV + iP * Dimension::NFF;
    //     bool* pf_cross = this->pf_cross + iP * Dimension::F;

    //     /////////////////////////////////////////////////
    //     for (int i = 0; i < Dimension::F; ++i) {
    //         double cross = 0.0e0;
    //         for (int j = 0; j < Dimension::N; ++j) {
    //             cross += p[j] * (grad[j] + dV[j * Dimension::FF + i * Dimension::Fadd1]);
    //         }
    //         pf_cross[i] = (cross > 0);
    //     }

    //     // bool group1 = false;
    //     // bool group2 = false;
    //     // for (int i = 0; i < Dimension::F; ++i) {
    //     //     if (pf_cross_last[i] && (!pf_cross[i])) {
    //     //         group1 = true;
    //     //     } else {
    //     //         group2 = true;
    //     //     }
    //     // }
    //     // if (group1 && group2) {  // do branching
    //     //     //
    //     // }
    // }
    return 0;
}

};  // namespace PROJECT_NS
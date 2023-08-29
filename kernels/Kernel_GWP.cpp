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

int Kernel_GWP::calc_Snuc(num_complex* Snuc,  // [P,P]
                          num_real* x,        // [P,N]
                          num_real* p,        // [P,N]
                          num_real* alpha,    // [N]
                          int P, int N) {
    for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
            num_complex term = 0.0e0;
            for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
                term += -0.25 * alpha[j] * (x[aj] - x[bj]) * (x[aj] - x[bj])   //
                        - 0.25 / alpha[j] * (p[aj] - p[bj]) * (p[aj] - p[bj])  //
                        + 0.5 * phys::math::im * (p[aj] + p[bj]) * (x[aj] - x[bj]);
            }
            Snuc[ab] = std::exp(term);
        }
    }
    return 0;
}

int Kernel_GWP::calc_Sele(num_complex* Sele,  // [P,P]
                          num_complex* c,     // [P,F]
                          int P, int F) {
    Eigen::Map<EigMX<num_complex>> Map_S(Sele, P, P);
    Eigen::Map<EigMX<num_complex>> Map_c(c, P, F);
    Map_S = Map_c.conjugate() * Map_c.transpose();
    return 0;
}

int Kernel_GWP::calc_dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                              num_real* x,            // [P,N]
                              num_real* p,            // [P,N]
                              num_real* m,            // [P,N]
                              num_real* f,            // [P,N]
                              num_real* alpha,        // [N]
                              int P, int N) {
    for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
            num_complex term = 0.0e0;
            for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
                term +=
                    p[bj] / m[bj] * 0.5 * alpha[j] * (x[aj] - x[bj] + (p[aj] + p[bj]) / (phys::math::im * alpha[j]));
                term += (-f[bj]) * 0.5 / alpha[j] * ((p[aj] - p[bj]) + phys::math::im * alpha[j] * (x[aj] - x[bj]));
                term += phys::math::im * p[bj] * p[bj] / (2 * m[bj]);
            }
            dtlnSnuc[ab] = term;
        }
    }
    return 0;
}

int Kernel_GWP::calc_dtSele(num_complex* dtSele,  // [P,P]
                            num_complex* c,       // [P,F]
                            num_complex* H,       // [P,F,F]
                            num_real* vpes,       // [P]
                            int P, int F) {
    int FF = F * F;
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

double Kernel_GWP::calc_density(num_complex* rhored,  // [F,F]
                                num_complex* Acoeff,  // [P]
                                num_complex* Snuc,    // [P,P]
                                num_complex* c,       // [P,F]
                                num_real xi,          //
                                num_real gamma,       //
                                int P, int F) {       //@bug
    Eigen::Map<EigMX<num_complex>> Map_rhored(rhored, Dimension::F, Dimension::F);
    Eigen::Map<EigMX<num_complex>> Map_Acoeff(Acoeff, Dimension::P, 1);
    Eigen::Map<EigMX<num_complex>> Map_Snuc(Snuc, Dimension::P, Dimension::P);
    Eigen::Map<EigMX<num_complex>> Map_c(c, Dimension::P, Dimension::F);
    Map_rhored =
        (xi * Map_c.adjoint() * Map_Acoeff.adjoint().asDiagonal() * Map_Snuc * Map_Acoeff.asDiagonal() * Map_c -
         gamma * (Map_Acoeff.adjoint() * Map_Snuc * Map_Acoeff).sum() * EigMXc::Identity(P, P)) /
        double(P);
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
    dt     = PM->get<double>("dt", LOC(), phys::time_d);  //
    alpha0 = PM->get<double>("alpha0", LOC(), 1.0f);      //
    gamma  = PM->get<double>("gamma", LOC(), 0.0f);       //
    xi     = 1 + Dimension::F * gamma;
}

void Kernel_GWP::init_data_impl(DataSet* DS) {
    alpha    = DS->reg<num_real>("integrator.alpha", Dimension::N);
    Acoeff   = DS->reg<num_complex>("integrator.Acoeff", Dimension::P);
    dtAcoeff = DS->reg<num_complex>("integrator.dtAcoeff", Dimension::P);
    Hcoeff   = DS->reg<num_complex>("integrator.Hcoeff", Dimension::PP);
    Hbasis   = DS->reg<num_complex>("integrator.Hbasis", Dimension::PP);
    rhored   = DS->reg<num_complex>("integrator.rhored", Dimension::FF);

    Snuc     = DS->reg<num_complex>("integrator.Snuc", Dimension::PP);
    Sele     = DS->reg<num_complex>("integrator.Sele", Dimension::PP);
    S        = DS->reg<num_complex>("integrator.S", Dimension::PP);
    invS     = DS->reg<num_complex>("integrator.invS", Dimension::PP);
    dtlnSnuc = DS->reg<num_complex>("integrator.dtlnSnuc", Dimension::PP);
    dtSele   = DS->reg<num_complex>("integrator.dtSele", Dimension::PP);
    //
    Udt  = DS->reg<num_complex>("integrator.Udt", Dimension::PFF);
    H    = DS->reg<num_complex>("model.rep.H", Dimension::PFF);
    vpes = DS->reg<num_real>("model.vpes", Dimension::P);
    grad = DS->reg<num_real>("model.grad", Dimension::PN);
    V    = DS->reg<num_real>("model.V", Dimension::PFF);
    dV   = DS->reg<num_real>("model.dV", Dimension::PNFF);
    E    = DS->reg<num_real>("model.rep.E", Dimension::PF);
    dE   = DS->reg<num_real>("model.rep.dE", Dimension::PNFF);
    x    = DS->reg<num_real>("integrator.x", Dimension::PN);
    p    = DS->reg<num_real>("integrator.p", Dimension::PN);
    m    = DS->reg<num_real>("integrator.m", Dimension::PN);
    f    = DS->reg<num_real>("integrator.f", Dimension::PN);
    c    = DS->reg<num_complex>("integrator.c", Dimension::PF);
}

void Kernel_GWP::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* w = Kernel_Elec::w + iP;
        num_complex* c = Kernel_Elec::c + iP * Dimension::F;
        num_complex* U = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc   = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0]     = 1.0e0;              ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;  ///< initial occupation
        /// < initial c (not used)
        for (int i = 0; i < Dimension::F; ++i) {
            double randu;
            Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
            c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
        }
        ARRAY_EYE(U, Dimension::F);  ///< initial propagator
    }
    for (int i = 0; i < Dimension::N; ++i) alpha[i] = alpha0;
    for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = 1.0e0;
    Kernel_Elec::c_init = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::T_init = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);
    exec_kernel(-1);
}

int Kernel_GWP::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* c_init  = Kernel_Elec::c_init + iP * Dimension::F;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;

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

        Kernel_Elec::ker_from_c(rho_nuc, c, xi, gamma, Dimension::F);
    }

    // calculate overlap integral & time derivative of overlap integral
    calc_Snuc(Snuc, x, p, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, Dimension::P, Dimension::F);
    calc_dtlnSnuc(dtlnSnuc, x, p, m, f, alpha, Dimension::P, Dimension::N);
    calc_dtSele(dtSele, c, H, vpes, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) S[ab] = Snuc[ab] * Sele[ab];
    // calculate inverse of overlap integral
    calc_invS(invS, S, Dimension::P);

    // ARRAY_SHOW(Snuc, Dimension::P, Dimension::P);
    // ARRAY_SHOW(Sele, Dimension::P, Dimension::P);
    // ARRAY_SHOW(dtlnSnuc, Dimension::P, Dimension::P);
    // ARRAY_SHOW(dtSele, Dimension::P, Dimension::P);
    // ARRAY_SHOW(S, Dimension::P, Dimension::P);
    // ARRAY_SHOW(invS, Dimension::P, Dimension::P);

    // calculate Hamiltonian between two configurations basis
    calc_Hbasis(Hbasis, vpes, grad, V, dV, x, p, m, alpha, Sele, c, Dimension::P, Dimension::N, Dimension::F);

    // ARRAY_SHOW(Hbasis, Dimension::P, Dimension::P);

    // calculate effective Hamiltonian
    for (int ab = 0; ab < Dimension::PP; ++ab) {
        Hcoeff[ab] = (Snuc[ab] * Hbasis[ab]                    //
                      - phys::math::im * S[ab] * dtlnSnuc[ab]  //
                      - phys::math::im * Snuc[ab] * dtSele[ab]);
    }

    ARRAY_MATMUL(Hcoeff, invS, Hcoeff, Dimension::P, Dimension::P, Dimension::P);

    // ARRAY_SHOW(Hcoeff, Dimension::P, Dimension::P);

    // update Acoeff
    ARRAY_MATMUL(dtAcoeff, Hcoeff, Acoeff, Dimension::P, Dimension::P, 1);
    if (stat >= 0)
        for (int a = 0; a < Dimension::P; ++a) { Acoeff[a] += -phys::math::im * dtAcoeff[a] * dt; }

    // ARRAY_SHOW(Acoeff, 1, Dimension::P);

    // calculate reduced density
    // ARRAY_SHOW(rhored, Dimension::F, Dimension::F);

    num_real scale = calc_density(rhored, Acoeff, Snuc, c, xi, gamma, Dimension::P, Dimension::F);

    // std::cout << scale << "\n";

    for (int ik = 0; ik < Dimension::FF; ++ik) rhored[ik] /= scale;
    scale = sqrt(scale);
    for (int a = 0; a < Dimension::P; ++a) Acoeff[a] /= scale;
    return 0;
}
};  // namespace PROJECT_NS
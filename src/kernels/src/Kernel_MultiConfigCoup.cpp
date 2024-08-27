#include "kids/Kernel_MultiConfigCoup.h"

#include <algorithm>

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_MultiConfigCoup::getName() { return "Kernel_MultiConfigCoup"; }

int Kernel_MultiConfigCoup::getType() const { return utils::hash(FUNCTION_NAME); }

int Kernel_MultiConfigCoup::calc_Ekin(span<kids_real> Ekin,  // [P]
                                      span<kids_real> p,     // [P,N]
                                      span<kids_real> m,     // [P,N]
                                      int P, int N) {
    for (int iP = 0; iP < P; ++iP) {
        span<kids_real> Ekin = Ekin.subspan(iP, 1);
        span<kids_real> p    = p.subspan(iP * N, N);
        span<kids_real> m    = m.subspan(iP * N, N);
        Ekin[0]              = 0;
        for (int j = 0; j < N; ++j) Ekin[0] += p[j] * p[j] / m[j];
        Ekin[0] /= 2;
    }
    return 0;
}

/**
 * the expression is exp(-0.25*a*(x1-x2)^2 -0.25*(p1-p2)/a + 0.5i*(p1+p2)(x1-x2) - i(g1-g2))
 */
int Kernel_MultiConfigCoup::calc_Snuc(span<kids_complex> Snuc,   // [P,P]
                                      span<kids_real>    x1,     // [P,N]
                                      span<kids_real>    p1,     // [P,N]
                                      span<kids_real>    m1,     // [P,N]
                                      span<kids_real>    g1,     // [P]
                                      span<kids_real>    x2,     // [P,N]
                                      span<kids_real>    p2,     // [P,N]
                                      span<kids_real>    m2,     // [P,N]
                                      span<kids_real>    g2,     // [P]
                                      span<kids_real>    alpha,  // [N]
                                      int P, int N) {
    // PRINT_ARRAY(Snuc, P, P);
    for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
            kids_complex term = 0.0e0;
            for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
                term += -0.25 * alpha[j] * (x1[aj] - x2[bj]) * (x1[aj] - x2[bj])         //
                        - 0.25 / alpha[j] * (p1[aj] - p2[bj]) * (p1[aj] - p2[bj])        //
                        + 0.5 * phys::math::im * (p1[aj] + p2[bj]) * (x1[aj] - x2[bj]);  //
            }
            Snuc[a * Dimension::P + b] = std::exp(term - phys::math::im * (g1[a] - g2[b]));
        }
    }
    return 0;
}

int Kernel_MultiConfigCoup::calc_Sele(span<kids_complex> Sele,   // [P,P]
                                      span<kids_complex> c1,     // [P,F]
                                      span<kids_complex> c2,     // [P,F]
                                      kids_real          xi,     //
                                      kids_real          gamma,  //
                                      int P, int F) {
    // Map_S = xi * Map_c1.conjugate() * Map_c2.transpose() - F * gamma * EigMXc::Identity(P, P);
    // @BAD!!! ? gamma should be zero
    for (int a = 0, ab = 0; a < Dimension::P; ++a) {
        span<kids_complex> c1a = c1.subspan(a * Dimension::F, Dimension::F);
        for (int b = 0; b < Dimension::P; ++b, ++ab) {
            span<kids_complex> c2b     = c2.subspan(b * Dimension::F, Dimension::F);
            Sele[a * Dimension::P + b] = xi * ARRAY_INNER_TRANS1(c1a.data(), c2b.data(), Dimension::P) -
                                         ((a == b) ? gamma : 0.0e0);  // ? F*gamma?
        }
    }
    return 0;  // @bug FATAL
}

int Kernel_MultiConfigCoup::calc_dtlnSnuc(span<kids_complex> dtlnSnuc,  // [P,P]
                                          span<kids_real>    x,         // [P,N]
                                          span<kids_real>    p,         // [P,N]
                                          span<kids_real>    m,         // [P,N]
                                          span<kids_real>    f,         // [P,N]
                                          span<kids_real>    alpha,     // [N]
                                          span<kids_real>    Ekin,      // [P]
                                          int P, int N) {
    // PRINT_ARRAY(dtlnSnuc, P, P);
    for (int a = 0, aN = 0, ab = 0; a < P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < P; ++b, ++ab, bN += N) {
            kids_complex term = 0.0e0;
            for (int j = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj) {
                term +=
                    p[bj] / m[bj] * 0.5 * alpha[j] * (x[aj] - x[bj] + (p[aj] + p[bj]) / (phys::math::im * alpha[j]));
                term += (-f[bj]) * 0.5 / alpha[j] * ((p[aj] - p[bj]) + phys::math::im * alpha[j] * (x[aj] - x[bj]));
                term += phys::math::im * p[bj] * p[bj] / (2 * m[bj]);
            }
            dtlnSnuc[a * Dimension::P + b] = term;
        }
    }
    // PRINT_ARRAY(dtlnSnuc, P, P);
    return 0;
}

int Kernel_MultiConfigCoup::calc_dtSele(span<kids_complex> dtSele,  // [P,P]
                                        span<kids_complex> Sele,    // [P,P]
                                        span<kids_complex> c,       // [P,F]
                                        span<kids_complex> H,       // [P,F,F]
                                        span<kids_real>    vpes,    // [P]
                                        int P, int F) {
    for (int a = 0, ab = 0, aF = 0; a < Dimension::P; ++a, aF += Dimension::F) {
        span<kids_complex> ca = c.subspan(aF, F);
        for (int b = 0, bF = 0, bFF = 0; b < Dimension::P; ++b, ++ab, bF += Dimension::F, bFF += Dimension::FF) {
            span<kids_complex> cb = c.subspan(bF, F);
            span<kids_complex> Hb = H.subspan(bFF, F * F);
            kids_complex term1    = ARRAY_INNER_VMV_TRANS1(ca.data(), H.data(), cb.data(), Dimension::F, Dimension::F);
            kids_complex term2    = ARRAY_INNER_TRANS1(ca.data(), cb.data(), Dimension::F);
            dtSele[a * Dimension::P + b] = -phys::math::im * (term1 + vpes[b] * term2);
        }
    }
    return 0;
}

double Kernel_MultiConfigCoup::calc_density(span<kids_complex> rhored,   // [F,F]
                                            span<kids_complex> Acoeff,   // [P]
                                            span<kids_complex> Snuc,     // [P,P]
                                            span<kids_complex> c,        // [P,F]
                                            span<kids_complex> Mtmp,     // [P,P]
                                            kids_real          xi,       // for kernel
                                            kids_real          gamma,    // for kernel
                                            int P_used, int P, int F) {  //@bug

    for (int a = 0, ab = 0; a < P; ++a) {
        for (int b = 0; b < P; ++b, ++ab) { Mtmp[ab] = std::conj(Acoeff[a]) * Acoeff[b] * Snuc[ab]; }
    }
    kids_complex val = ARRAY_INNER_VMV_TRANS1(Acoeff.data(), Snuc.data(), Acoeff.data(), P, P);
    ARRAY_MATMUL3_TRANS1(rhored.data(), c.data(), Mtmp.data(), c.data(), F, P, P, F);
    kids_real trace = 0.0e0;
    for (int i = 0, ik = 0; i < F; ++i) {
        for (int k = 0; k < F; ++k, ++ik) {
            rhored[ik] = xi * rhored[ik];
            if (i == k) {
                rhored[ik] -= gamma * val;
                trace += std::real(rhored[ik]);
            }
        }
    }
    return trace;
}

int Kernel_MultiConfigCoup::calc_Hbasis(span<kids_complex> Hbasis,  // [P,P]
                                        span<kids_real>    vpes,    // [P]
                                        span<kids_real>    grad,    // [P,N]
                                        span<kids_real>    V,       // [P,F,F]
                                        span<kids_real>    dV,      // [P,N,F,F]
                                        span<kids_real>    x,       // [P,N]
                                        span<kids_real>    p,       // [P,N]
                                        span<kids_real>    m,       // [P,N]
                                        span<kids_real>    alpha,   // [N]
                                        span<kids_complex> Sele,    // [P,P]
                                        span<kids_complex> c,       // [P,F]
                                        int P, int N, int F) {
    int FF  = F * F;
    int NFF = N * FF;
    for (int a = 0, aN = 0, aF = 0, ab = 0; a < Dimension::P; ++a, aN += N, aF += F) {
        span<kids_real>    Va = V.subspan(a * FF, FF);
        span<kids_complex> ca = c.subspan(aF, F);
        for (int b = 0, bN = 0, bF = 0; b < Dimension::P; ++b, ++ab, bN += N, bF += F) {
            span<kids_real>    Vb = V.subspan(b * FF, FF);
            span<kids_complex> cb = c.subspan(bF, F);

            kids_complex Tab = 0.0e0;
            kids_complex Vab = 0.0e0;

            Vab += 0.5e0 * (vpes[a] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca.data(), Va.data(), cb.data(), F, F));
            Vab += 0.5e0 * (vpes[b] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca.data(), Vb.data(), cb.data(), F, F));
            for (int j = 0, jFF = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj, jFF += FF) {
                span<kids_real> dVaj = dV.subspan(aj * FF, FF);
                span<kids_real> dVbj = dV.subspan(bj * FF, FF);
                kids_complex    xabj = 0.5e0 * (x[aj] + x[bj] + (p[aj] - p[bj]) / (phys::math::im * alpha[j]));
                kids_complex    pabj = 0.5e0 * (p[aj] + p[bj] + (phys::math::im * alpha[j]) * (x[aj] - x[bj]));
                Tab += alpha[j] / (4 * m[bj]) + pabj * pabj / (2 * m[bj]);

                Vab += 0.5e0 * (xabj - x[aj]) *
                       (grad[aj] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca.data(), dVaj.data(), cb.data(), F, F));
                Vab += 0.5e0 * (xabj - x[bj]) *
                       (grad[bj] * Sele[ab] + ARRAY_INNER_VMV_TRANS1(ca.data(), dVbj.data(), cb.data(), F, F));
            }
            Hbasis[ab] = Tab + Vab;
        }
    }
    return 0;
};

int Kernel_MultiConfigCoup::calc_Hbasis_adia(span<kids_complex> Hbasis,  // [P,P]
                                             span<kids_real>    E,       // [P,F]
                                             span<kids_real>    dE,      // [P,N,F,F]
                                             span<kids_real>    x,       // [P,N]
                                             span<kids_real>    p,       // [P,N]
                                             span<kids_real>    m,       // [P,N]
                                             span<kids_real>    alpha,   // [N]
                                             span<kids_complex> c,       // [P,F]
                                             int P, int N, int F) {
    int FF  = F * F;
    int NFF = N * FF;
    for (int a = 0, aN = 0, ab = 0; a < Dimension::P; ++a, aN += N) {
        for (int b = 0, bN = 0; b < Dimension::P; ++b, ++ab, bN += N) {
            span<kids_real>    Ea  = E.subspan(a * F, F);
            span<kids_real>    Eb  = E.subspan(b * F, F);
            span<kids_real>    dEa = dE.subspan(a * NFF, NFF);
            span<kids_real>    dEb = dE.subspan(b * NFF, NFF);
            span<kids_complex> ca  = c.subspan(a * F, F);
            span<kids_complex> cb  = c.subspan(b * F, F);

            kids_complex Tab = 0.0e0;
            kids_complex Eab = 0.0e0;

            for (int j = 0, jik = 0, jFF = 0, aj = aN, bj = bN; j < N; ++j, ++aj, ++bj, jFF += Dimension::FF) {
                span<kids_real> dEaj = dEa.subspan(jFF, FF);
                span<kids_real> dEbj = dEb.subspan(jFF, FF);
                kids_complex    xabj = 0.5e0 * (x[aj] + x[bj] + (p[aj] - p[bj]) / (phys::math::im * alpha[j]));
                kids_complex    pabj = 0.5e0 * (p[aj] + p[bj] + (phys::math::im * alpha[j]) * (x[aj] - x[bj]));
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

void Kernel_MultiConfigCoup::setInputParam_impl(std::shared_ptr<Param> PM) {
    dt            = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d);  //
    alpha0        = _param->get_real({"solver.alpha0"}, LOC(), 1.0f);                  //
    width_scaling = _param->get_real({"solver.width_scaling"}, LOC(), 1.0f);           //
    break_thres   = _param->get_real({"solver.break_thres"}, LOC(), 1.0f);             //
    P_used0       = _param->get_int({"solver.P_initial"}, LOC(), 1);                   //
    max_clone     = _param->get_int({"solver.max_clone"}, LOC(), 0);                   //
    gamma         = _param->get_real({"solver.gamma"}, LOC(), 0.0f);                   //
    xi            = 1 + Dimension::F * gamma;

    impl_type          = _param->get_int({"solver.impl_type"}, LOC(), 0);           //
    samp_type          = _param->get_int({"solver.samp_type"}, LOC(), 0);           //
    aset_type          = _param->get_int({"solver.aset_type"}, LOC(), 0);           //
    time_displace_step = _param->get_int({"solver.time_displace_step"}, LOC(), 1);  //
}

void Kernel_MultiConfigCoup::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    alpha    = DS->def(DATA::integrator::alpha);  //? N size check again
    Xcoeff   = DS->def(DATA::integrator::Xcoeff);
    Acoeff   = DS->def(DATA::integrator::Acoeff);
    dtAcoeff = DS->def(DATA::integrator::dtAcoeff);
    Hcoeff   = DS->def(DATA::integrator::Hcoeff);
    Hbasis   = DS->def(DATA::integrator::Hbasis);
    UXdt     = DS->def(DATA::integrator::UXdt);
    UYdt     = DS->def(DATA::integrator::UYdt);
    rhored   = DS->def(DATA::integrator::rhored);
    rhored2  = DS->def(DATA::integrator::rhored2);
    rhored3  = DS->def(DATA::integrator::rhored3);

    Snuc     = DS->def(DATA::integrator::Snuc);
    Sele     = DS->def(DATA::integrator::Sele);
    S        = DS->def(DATA::integrator::S);
    invS     = DS->def(DATA::integrator::invS);
    dtlnSnuc = DS->def(DATA::integrator::dtlnSnuc);
    dtSele   = DS->def(DATA::integrator::dtSele);
    L        = DS->def(DATA::integrator::GWP::L);
    L1       = DS->def(DATA::integrator::GWP::L1);
    L2       = DS->def(DATA::integrator::GWP::L2);
    R        = DS->def(DATA::integrator::GWP::R);
    R1       = DS->def(DATA::integrator::GWP::R1);
    R2       = DS->def(DATA::integrator::GWP::R2);
    S1       = DS->def(DATA::integrator::GWP::S1);
    S1h      = DS->def(DATA::integrator::GWP::S1h);
    invS1h   = DS->def(DATA::integrator::GWP::invS1h);
    S2       = DS->def(DATA::integrator::GWP::S2);
    S2h      = DS->def(DATA::integrator::GWP::S2h);
    invS2h   = DS->def(DATA::integrator::GWP::invS2h);
    Sx       = DS->def(DATA::integrator::GWP::Sx);

    Ekin          = DS->def(DATA::integrator::Ekin);
    g             = DS->def(DATA::integrator::g);
    clone_account = DS->def(DATA::integrator::clone_account);
    // pf_cross      = DS->def<kids_bool>(DATA::integrator::pf_cross);

    //
    Udt  = DS->def(DATA::integrator::Udt);
    H    = DS->def(DATA::model::rep::H);
    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    eig  = DS->def(DATA::model::rep::eig);  // check change to eig
    T    = DS->def(DATA::model::rep::T);
    dE   = DS->def(DATA::model::rep::dE);
    x    = DS->def(DATA::integrator::x);
    p    = DS->def(DATA::integrator::p);
    m    = DS->def(DATA::integrator::m);
    f    = DS->def(DATA::integrator::f);
    c    = DS->def(DATA::integrator::c);

    w       = DS->def(DATA::integrator::w);
    U       = DS->def(DATA::integrator::U);
    occ_nuc = DS->def(DATA::integrator::occ_nuc);
    c_init  = DS->def(DATA::init::c);
    T_init  = DS->def(DATA::init::T);

    fun_diag_F = DS->def(DATA::integrator::tmp::fun_diag_F);
    fun_diag_P = DS->def(DATA::integrator::tmp::fun_diag_P);
    MatR_PP    = DS->def(DATA::integrator::tmp::MatR_PP);
    MatC_PP    = DS->def(DATA::integrator::tmp::MatC_PP);
    I_PP       = DS->def(DATA::integrator::tmp::I_PP);
    Ubranch    = DS->def(DATA::integrator::Ubranch);

    x_last    = DS->def(DATA::last::x);
    p_last    = DS->def(DATA::last::p);
    grad_last = DS->def(DATA::last::grad);
    dV_last   = DS->def(DATA::last::dV);
    g_last    = DS->def(DATA::last::g);
    c_last    = DS->def(DATA::last::c);

    P_used_ptr = DS->def(DATA::integrator::P_used);
    norm_ptr   = DS->def(DATA::integrator::norm);
    veF        = DS->def(DATA::integrator::veF);
    ve         = DS->def(DATA::integrator::ve);
}

Status& Kernel_MultiConfigCoup::initializeKernel_impl(Status& stat) { return stat; }

// {
//     // @begin debug
//     // for (int iP = 0; iP < Dimension::P; ++iP) {
//     //     span<kids_real> x = this->x + iP * Dimension::N;
//     //     span<kids_real> p = this->p + iP * Dimension::N;
//     //     for (int j = 0; j < Dimension::N; ++j) {
//     //         x[j] = iP * 0.02 * (iP % 2 - 0.5) + 0.1 * j;
//     //         p[j] = -iP * 0.02 * (iP % 2 - 0.5) + 0.2 * j;
//     //     }
//     // }
//     // @end debug

//     if (samp_type < 3) {  // overlap or neighbourhood re-sampling
//         for (int iP = 0; iP < Dimension::P; ++iP) {
//             span<kids_complex> w       = this->w + iP;
//             span<kids_complex> c       = this->c + iP * Dimension::F;
//             span<kids_complex> U       = this->U + iP * Dimension::FF;
//             int*          occ_nuc = this->occ_nuc + iP;

//             span<kids_real> x = this->x + iP * Dimension::N;
//             span<kids_real> p = this->p + iP * Dimension::N;

//             /////////////////////////////////////////////////////////////////
//             if (samp_type == 1)
//                 for (int j = 0; j < Dimension::N; ++j) {
//                     x[j] = this->x[j];
//                     p[j] = this->p[j];
//                 }

//             if (samp_type == 2 && iP > 0) {
//                 for (int j = 0; j < Dimension::N; ++j) {
//                     double randu;
//                     Kernel_Random::rand_gaussian(&randu);
//                     randu = sqrt(iP * iP % (j + 1));
//                     x[j]  = this->x[j] + width_scaling * randu / sqrt(Dimension::N * alpha0);
//                     Kernel_Random::rand_gaussian(&randu);
//                     randu = sqrt(iP % (j + 1));
//                     p[j]  = this->p[j] + width_scaling * randu * sqrt(alpha0 / Dimension::N);
//                 }
//             }

//             *w       = 1.0e0;  ///< initial measure
//             *occ_nuc = occ0;   ///< initial occupation
//             /// < initial c (not used)
//             for (int i = 0; i < Dimension::F; ++i) {
//                 double randu;
//                 Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
//                 c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
//             }
//             ARRAY_EYE(U, Dimension::F);  ///< initial propagator
//         }
//     }
//     if (samp_type == 3) {  /// time displaced re-sampling
//         w[0]       = 1.0e0;
//         occ_nuc[0] = occ0;
//         for (int i = 0; i < Dimension::F; ++i) {
//             double randu;
//             Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
//             c[i] = (i == occ_nuc[0] ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
//         }
//         ARRAY_EYE(U, Dimension::F);
//         for (int iP = 1; iP < Dimension::P; ++iP) {
//             span<kids_real>    x_now       = x + iP * Dimension::N;
//             span<kids_real>    p_now       = p + iP * Dimension::N;
//             span<kids_real>    f_now       = f + iP * Dimension::N;
//             span<kids_complex> U_now       = this->U + iP * Dimension::FF;
//             span<kids_complex> c_now       = this->c + iP * Dimension::F;
//             span<kids_complex> rho_nuc_now = this->rho_nuc + iP * Dimension::FF;

//             span<kids_real>    x_prev = x + std::max({iP - 2, 0}) * Dimension::N;
//             span<kids_real>    p_prev = p + std::max({iP - 2, 0}) * Dimension::N;
//             span<kids_real>    f_prev = f + std::max({iP - 2, 0}) * Dimension::N;
//             span<kids_complex> U_prev = this->U + std::max({iP - 2, 0}) * Dimension::FF;

//             span<kids_real>    eig_now = eig + iP * Dimension::F;
//             span<kids_real>    T_now   = T + iP * Dimension::FF;
//             span<kids_complex> Udt_now = Udt + iP * Dimension::FF;

//             kids_real signdt = (iP % 2 == 0) ? dt : -dt;

//             for (int j = 0; j < Dimension::N; ++j) x_now[j] = x_prev[j], p_now[j] = p_prev[j], f_now[j] = f_prev[j];

//             for (int istep_displace = 0; istep_displace < time_displace_step; ++istep_displace) {
//                 for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
//                 for (int j = 0; j < Dimension::N; ++j) x_now[j] += p_now[j] / m[j] * signdt;
//                 _kmodel->executeKernel(stat);  // only iP needed
//                 _krepr->executeKernel(stat);   // only iP needed
//                 switch (Kernel_Representation::ele_repr_type) {
//                     case RepresentationPolicy::Diabatic: {
//                         for (int i = 0; i < Dimension::F; ++i)
//                             fun_diag_F[i] = exp(-phys::math::im * eig_now[i] * signdt);
//                         ARRAY_MATMUL3_TRANS2(Udt_now, T_now, fun_diag_F, T_now, Dimension::F, Dimension::F, 0,
//                                              Dimension::F);
//                         break;
//                     }
//                 }
//                 ARRAY_MATMUL(U_now, Udt_now, U_prev, Dimension::F, Dimension::F, Dimension::F);
//                 ARRAY_MATMUL(c_now, U_now, c, Dimension::F, Dimension::F, 1);
//                 elec_utils::ker_from_c(rho_nuc_now, c_now, xi, gamma, Dimension::F);
//                 _kforce->executeKernel(stat);
//                 for (int j = 0; j < Dimension::N; ++j) p_now[j] -= f_now[j] * 0.5 * signdt;
//             }
//             // PRINT_ARRAY(x_now, 1, Dimension::N);
//             // PRINT_ARRAY(p_now, 1, Dimension::N);
//         }

//         for (int iP = 0; iP < Dimension::P; ++iP) {
//             span<kids_complex> w       = this->w + iP;
//             span<kids_complex> c       = this->c + iP * Dimension::F;
//             span<kids_complex> U       = this->U + iP * Dimension::FF;
//             int*          occ_nuc = this->occ_nuc + iP;

//             /////////////////////////////////////////////////////////////////

//             w[0]       = 1.0e0;  ///< initial measure
//             occ_nuc[0] = occ0;   ///< initial occupation
//             /// < use time-displaced c ?
//             // for (int i = 0; i < Dimension::F; ++i) {
//             //     double randu;
//             //     Kernel_Random::rand_uniform(&randu, 1, phys::math::twopi);
//             //     c[i] = (i == *occ_nuc ? sqrt(1 + gamma) : sqrt(gamma)) * exp(phys::math::im * randu) / sqrt(xi);
//             // }
//             ARRAY_EYE(U, Dimension::F);  ///< initial propagator reset to identity
//         }
//     }

//     _kmodel->executeKernel(stat);
//     _krepr->executeKernel(stat);

//     for (int i = 0; i < Dimension::N; ++i) alpha[i] = alpha0;
//     for (int a = 0; a < Dimension::P; ++a) g[a] = 0.0e0;
//     if (aset_type == 0) {
//         for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = (a == 0) ? 1.0e0 : 0.0e0;
//     } else if (aset_type == 1) {
//         for (int a = 0; a < Dimension::P; ++a) Acoeff[a] = 1.0e0;  // @only for test
//     }

//     // normalization of A
//     P_used        = P_used0;
//     P_used_ptr[0] = P_used;
//     // std::cout << "P_used_INIT = " << P_used << "\n";
//     for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
//     calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
//     calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
//     for (int ab = 0; ab < Dimension::PP; ++ab) { S[ab] = Snuc[ab] * Sele[ab]; }
//     kids_complex scale;
//     ARRAY_MATMUL3_TRANS1(&scale, Acoeff, S, Acoeff, 1, Dimension::P, Dimension::P, 1);
//     for (int a = 0; a < Dimension::P; ++a) Acoeff[a] /= sqrt(abs(scale));
//     ARRAY_MATMUL3_TRANS1(&scale, Acoeff, Snuc, Acoeff, 1, Dimension::P, Dimension::P, 1);
//     xi = 1.0e0 + Dimension::F * gamma * std::abs(scale);

//     std::cout << "init scale: " << scale << "\n";
//     std::cout << "init xi: " << xi << "\n";

//     norm_ptr[0] = 1.0e0;

//     // PRINT_ARRAY(Acoeff, 1, Dimension::P);
//     for (int a = 0; a < Dimension::P; ++a) g_last[a] = g[a];
//     for (int aj = 0; aj < Dimension::PN; ++aj) x_last[aj] = x[aj];
//     for (int aj = 0; aj < Dimension::PN; ++aj) p_last[aj] = p[aj];
//     for (int ai = 0; ai < Dimension::PF; ++ai) c_last[ai] = c[ai];

//     double dt_backup = dt;
//     dt               = 0.0e0;
//     executeKernel(stat);  // don't evolve!!!
//     dt = dt_backup;

//     return stat;
// }

Status& Kernel_MultiConfigCoup::impl_0(Status& stat) {
    // calculate overlap integral & time derivative of overlap integral

    // PRINT_ARRAY(p, Dimension::P, Dimension::N);
    // PRINT_ARRAY(x, Dimension::P, Dimension::N);
    // PRINT_ARRAY(c, Dimension::P, Dimension::F);
    // PRINT_ARRAY(g, 1, Dimension::P);

    calc_Ekin(Ekin, p, m, Dimension::P, Dimension::N);
    calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
    calc_dtlnSnuc(dtlnSnuc, x, p, m, f, alpha, Ekin, Dimension::P, Dimension::N);
    calc_dtSele(dtSele, Sele, c, H, vpes, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) { S[ab] = Snuc[ab] * Sele[ab]; }
    // calculate inverse of overlap integral
    EigenSolve(L1.data(), R1.data(), S.data(), Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / L1[a];
    ARRAY_MATMUL3_TRANS2(invS.data(), R1.data(), fun_diag_P.data(), R1.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);

    // calculate Hamiltonian between two configurations basis
    calc_Hbasis(Hbasis, vpes, grad, V, dV, x, p, m, alpha, Sele, c, Dimension::P, Dimension::N, Dimension::F);
    // calculate effective Hamiltonian
    for (int ab = 0; ab < Dimension::PP; ++ab) {
        Hcoeff[ab] =
            (Snuc[ab] * Hbasis[ab] - phys::math::im * S[ab] * dtlnSnuc[ab] - phys::math::im * Snuc[ab] * dtSele[ab]);
    }
    ARRAY_MATMUL(Hcoeff.data(), invS.data(), Hcoeff.data(), Dimension::P, Dimension::P, Dimension::P);

    // PRINT_ARRAY(Hcoeff, Dimension::P, Dimension::P);
    ARRAY_EXP_MAT_GENERAL(UXdt.data(), Hcoeff.data(), -phys::math::im * dt, Dimension::P);

    // PRINT_ARRAY(UXdt, Dimension::P, Dimension::P);

    // update Acoeff
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    ARRAY_MATMUL(Acoeff.data(), UXdt.data(), Acoeff.data(), Dimension::P, Dimension::P, 1);
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    kids_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff.data(), S.data(), Acoeff.data(), 1, Dimension::P, Dimension::P, 1);
    norm_ptr[0] *= std::abs(scale);
    for (int a = 0; a < Dimension::P; ++a) Acoeff[a] /= sqrt(abs(scale));

    for (int a = 0; a < Dimension::P; ++a) g[a] += Ekin[a] * dt;

    cloning();
    death();
    return stat;
}

Status& Kernel_MultiConfigCoup::impl_1(Status& stat) {
    // calculate overlap integral
    calc_Ekin(Ekin, p, m, Dimension::P, Dimension::N);
    calc_Snuc(Snuc, x_last, p_last, m, g_last, x_last, p_last, m, g_last, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c_last, c_last, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) S1[ab] = Snuc[ab] * Sele[ab];

    // PRINT_ARRAY(Acoeff, 1, Dimension::P);
    // PRINT_ARRAY(S1, Dimension::P, Dimension::P);

    // PRINT_ARRAY(S1, Dimension::P, Dimension::P);

    calc_Snuc(Snuc, x, p, m, g, x_last, p_last, m, g_last, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c_last, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) Sx[ab] = Snuc[ab] * Sele[ab];

    // PRINT_ARRAY(Sx, Dimension::P, Dimension::P);

    calc_Snuc(Snuc, x, p, m, g, x, p, m, g, alpha, Dimension::P, Dimension::N);
    calc_Sele(Sele, c, c, 1, 0, Dimension::P, Dimension::F);
    for (int ab = 0; ab < Dimension::PP; ++ab) S2[ab] = Snuc[ab] * Sele[ab];
    for (int ab = 0; ab < Dimension::PP; ++ab) S[ab] = S2[ab];

    // PRINT_ARRAY(S2, Dimension::P, Dimension::P);

    EigenSolve(L1.data(), R1.data(), S1.data(), Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = sqrt(L1[a]);
    ARRAY_MATMUL3_TRANS2(S1h.data(), R1.data(), fun_diag_P.data(), R1.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / sqrt(L1[a]);
    ARRAY_MATMUL3_TRANS2(invS1h.data(), R1.data(), fun_diag_P.data(), R1.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);

    EigenSolve(L2.data(), R2.data(), S2.data(), Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = sqrt(L2[a]);
    ARRAY_MATMUL3_TRANS2(S2h.data(), R2.data(), fun_diag_P.data(), R2.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = 1.0e0 / sqrt(L2[a]);
    ARRAY_MATMUL3_TRANS2(invS2h.data(), R2.data(), fun_diag_P.data(), R2.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);

    // calculate Hamiltonian between two configurations basis
    calc_Hbasis(Hbasis, vpes, grad, V, dV, x, p, m, alpha, Sele, c, Dimension::P, Dimension::N, Dimension::F);
    ARRAY_MATMUL(Hcoeff.data(), invS1h.data(), Hbasis.data(), Dimension::P, Dimension::P, Dimension::P);
    ARRAY_MATMUL(Hcoeff.data(), Hcoeff.data(), invS1h.data(), Dimension::P, Dimension::P, Dimension::P);

    // PRINT_ARRAY(Hcoeff, Dimension::P, Dimension::P);

    EigenSolve(L.data(), R.data(), Hcoeff.data(), Dimension::P);
    for (int a = 0; a < Dimension::P; ++a) fun_diag_P[a] = exp(-phys::math::im * L[a] * dt);
    ARRAY_MATMUL3_TRANS2(UXdt.data(), R.data(), fun_diag_P.data(), R.data(), Dimension::P, Dimension::P, 0,
                         Dimension::P);

    ARRAY_MATMUL(UYdt.data(), Sx.data(), invS1h.data(), Dimension::P, Dimension::P, Dimension::P);
    ARRAY_MATMUL(UYdt.data(), invS2h.data(), UYdt.data(), Dimension::P, Dimension::P, Dimension::P);
    ARRAY_CORRECT_U(UYdt.data(), Dimension::P);

    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;
    ARRAY_MATMUL(Xcoeff.data(), S1h.data(), Acoeff.data(), Dimension::P, Dimension::P, 1);

    kids_complex cnorm;
    ARRAY_MATMUL_TRANS1(&cnorm, Xcoeff.data(), Xcoeff.data(), 1, Dimension::P, 1);
    std::cout << "norm 1 = " << cnorm << "\n";

    ARRAY_MATMUL(Xcoeff.data(), UXdt.data(), Xcoeff.data(), Dimension::P, Dimension::P, 1);

    ARRAY_MATMUL_TRANS1(&cnorm, Xcoeff.data(), Xcoeff.data(), 1, Dimension::P, 1);
    std::cout << "norm 2 = " << cnorm << "\n";

    ARRAY_MATMUL(Xcoeff.data(), UYdt.data(), Xcoeff.data(), Dimension::P, Dimension::P, 1);

    ARRAY_MATMUL_TRANS1(&cnorm, Xcoeff.data(), Xcoeff.data(), 1, Dimension::P, 1);
    std::cout << "norm 3 = " << cnorm << "\n";

    // ARRAY_MATMUL(Xcoeff, invS1h, Xcoeff, Dimension::P, Dimension::P, 1);
    // ARRAY_MATMUL(Xcoeff, Sx, Xcoeff, Dimension::P, Dimension::P, 1);
    // ARRAY_MATMUL(Xcoeff, invS2h, Xcoeff, Dimension::P, Dimension::P, 1);
    ARRAY_MATMUL(Acoeff.data(), invS2h.data(), Xcoeff.data(), Dimension::P, Dimension::P, 1);
    for (int a = P_used; a < Dimension::P; ++a) Acoeff[a] = 0.0e0;

    std::cout << "P_used = " << P_used << "\n";
    // PRINT_ARRAY(Acoeff, 1, Dimension::P);
    // PRINT_ARRAY(S2, Dimension::P, Dimension::P);
    cloning();
    death();

    for (int a = 0; a < Dimension::P; ++a) g_last[a] = g[a];
    for (int aj = 0; aj < Dimension::PN; ++aj) x_last[aj] = x[aj];
    for (int aj = 0; aj < Dimension::PN; ++aj) p_last[aj] = p[aj];
    for (int ai = 0; ai < Dimension::PF; ++ai) c_last[ai] = c[ai];
    // update phase
    for (int a = 0; a < Dimension::P; ++a) g[a] += Ekin[a] * dt;
    return stat;
}
Status& Kernel_MultiConfigCoup::executeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        span<kids_complex> U       = this->U.subspan(iP * Dimension::FF, Dimension::FF);
        span<kids_complex> c       = this->c.subspan(iP * Dimension::F, Dimension::F);
        span<kids_complex> c_init  = this->c_init.subspan(iP * Dimension::F, Dimension::F);
        span<kids_complex> rho_nuc = this->rho_nuc.subspan(iP * Dimension::FF, Dimension::FF);
        span<kids_real>    T       = this->T.subspan(iP * Dimension::FF, Dimension::FF);
        span<kids_real>    T_init  = this->T_init.subspan(iP * Dimension::FF, Dimension::FF);

        /////////////////////////////////////////////////////////////////

        for (int i = 0; i < Dimension::F; ++i) c[i] = c_init[i];
        // 1) transform from inp_repr => ele_repr
        Kernel_Representation::transform(c.data(), T_init.data(), Dimension::F,  //
                                         Kernel_Representation::inp_repr_type,   //
                                         Kernel_Representation::ele_repr_type,   //
                                         SpacePolicy::H);
        // 2) propagte along ele_repr
        ARRAY_MATMUL(c.data(), U.data(), c.data(), Dimension::F, Dimension::F, 1);
        // 3) transform back from ele_repr => inp_repr
        Kernel_Representation::transform(c.data(), T.data(), Dimension::F,      //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::H);
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
    kids_complex scale;
    ARRAY_MATMUL3_TRANS1(&scale, Acoeff.data(), Snuc.data(), Acoeff.data(), 1, Dimension::P, Dimension::P, 1);
    std::cout << "t scale : " << scale << "\n";
    xi = 1.0e0 + Dimension::F * gamma * std::abs(scale);
    calc_density(rhored, Acoeff, Snuc, c, MatC_PP, xi, gamma, P_used, Dimension::P, Dimension::F);

    ARRAY_CLEAR(rhored2.data(), Dimension::FF);
    for (int a = 0; a < P_used; ++a) {
        span<kids_complex> ca = c.subspan(a * Dimension::F, Dimension::F);
        for (int i = 0, ik = 0; i < Dimension::F; ++i) {
            for (int k = 0; k < Dimension::F; ++k, ++ik) {
                rhored2[ik] += xi * ca[i] * std::conj(ca[k]) - ((i == k) ? gamma : 0.0e0);
            }
        }
    }
    for (int ik = 0; ik < Dimension::FF; ++ik) rhored2[ik] /= double(P_used);
    return stat;
}

int Kernel_MultiConfigCoup::cloning() {
    P_used_ptr[0] = P_used;
    if (P_used >= Dimension::P) return 0;
    int P_increase = P_used;

    // calc_Snuc(Snuc, this->x, this->p, this->m, this->g, this->x, this->p, this->m, this->g, alpha, Dimension::P,
    //           Dimension::N);
    // calc_density(rhored, Acoeff, Snuc, c, xi, gamma, P_used, Dimension::P, Dimension::F);
    // PRINT_ARRAY(Snuc, Dimension::P, Dimension::P);
    // PRINT_ARRAY(c, Dimension::P, Dimension::F);
    // PRINT_ARRAY(Acoeff, Dimension::P, 1);
    // PRINT_ARRAY(rhored, Dimension::F, Dimension::F);

    // calc_density(rhored, Acoeff, Snuc, c, 1, 0, P_used, Dimension::P, Dimension::F);
    // PRINT_ARRAY(rhored, Dimension::F, Dimension::F);

    for (int iP = 0; iP < P_used; ++iP) {
        span<kids_real>    g      = this->g.subspan(iP, 1);
        span<kids_real>    x      = this->x.subspan(iP * Dimension::N, Dimension::N);
        span<kids_real>    p      = this->p.subspan(iP * Dimension::N, Dimension::N);
        span<kids_real>    f      = this->f.subspan(iP * Dimension::N, Dimension::N);
        span<kids_real>    grad   = this->grad.subspan(iP * Dimension::N, Dimension::N);
        span<kids_real>    V      = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        span<kids_real>    dV     = this->dV.subspan(iP * Dimension::NFF, Dimension::NFF);
        span<kids_complex> c      = this->c.subspan(iP * Dimension::F, Dimension::F);
        span<kids_complex> c_init = this->c_init.subspan(iP * Dimension::F, Dimension::F);
        span<kids_complex> U      = this->U.subspan(iP * Dimension::FF, Dimension::FF);

        /////////////////////////////////////////////////

        double max_break_val = 0.0e0;
        int    break_state_i = -1, break_state_k = -1;
        // Eigen::Map<EigMXr> Map_veF(veF, Dimension::FF, 1);
        // Eigen::Map<EigMXr> Map_dV(dV, Dimension::N, Dimension::FF);
        // Eigen::Map<EigMXr> Map_p(p, Dimension::N, 1);
        // Eigen::Map<EigMXr> Map_m(m, Dimension::N, 1);
        // Map_veF = Map_dV.transpose() * (Map_p.array() / Map_m.array()).matrix();
        for (int j = 0; j < Dimension::N; ++j) ve[j] = p[j] / m[j];
        ARRAY_MATMUL(veF.data(), ve.data(), dV.data(), 1, Dimension::N, Dimension::FF);


        for (int i = 0; i < Dimension::F; ++i) {
            for (int k = i + 1; k < Dimension::F; ++k) {
                double break_val =
                    std::abs(c[i] * c[i] * c[k] * c[k] * (veF[i * Dimension::Fadd1] - veF[k * Dimension::Fadd1]) /
                             (1e-6 + std::abs(V[i * Dimension::Fadd1] - V[k * Dimension::Fadd1])));
                if (break_val > max_break_val) {
                    max_break_val = break_val;
                    break_state_i = i;
                    break_state_k = k;
                }
            }
        }
        if (max_break_val > break_thres) {
            int break_state = (std::abs(c[break_state_i]) > std::abs(c[break_state_k])) ? break_state_i : break_state_k;
            double       norm_b       = std::abs(c[break_state]);
            double       norm_a       = sqrt(1 - norm_b * norm_b);
            kids_complex norm_phase_a = norm_a * (norm_a + phys::math::im * norm_b);
            kids_complex norm_phase_b = norm_b * (norm_b - phys::math::im * norm_a);

            // std::cout << "norm_a, norm_b, iP, P_used, P_increase: " << norm_a << ", " << norm_b << ", " << iP << ","
            //           << P_used << ", " << P_increase << "\n";

            span<kids_real>    g_new      = this->g.subspan(P_increase, 1);
            span<kids_real>    x_new      = this->x.subspan(P_increase * Dimension::N, Dimension::N);
            span<kids_real>    p_new      = this->p.subspan(P_increase * Dimension::N, Dimension::N);
            span<kids_real>    f_new      = this->f.subspan(P_increase * Dimension::N, Dimension::N);
            span<kids_real>    grad_new   = this->grad.subspan(P_increase * Dimension::N, Dimension::N);
            span<kids_real>    dV_new     = this->dV.subspan(P_increase * Dimension::NFF, Dimension::NFF);
            span<kids_complex> c_new      = this->c.subspan(P_increase * Dimension::F, Dimension::F);
            span<kids_complex> c_init_new = this->c_init.subspan(P_increase * Dimension::F, Dimension::F);
            span<kids_complex> U_new      = this->U.subspan(P_increase * Dimension::FF, Dimension::FF);

            g_new[0] = g[0];
            for (int j = 0; j < Dimension::N; ++j) x_new[j] = x[j];
            for (int j = 0; j < Dimension::N; ++j) p_new[j] = p[j];
            for (int j = 0; j < Dimension::N; ++j) f_new[j] = f[j];
            for (int j = 0; j < Dimension::N; ++j) grad_new[j] = grad[j];
            for (int j = 0; j < Dimension::NFF; ++j) dV_new[j] = dV[j];

            for (int i = 0; i < Dimension::F; ++i) fun_diag_F[i] = c[i];
            for (int i = 0; i < Dimension::F; ++i) c_new[i] = ((i == break_state) ? 0.0e0 : (c[i] / norm_phase_a));
            for (int i = 0; i < Dimension::F; ++i) c_init_new[i] = c_init[i];
            ARRAY_MATMUL_TRANS2(Ubranch.data(), c_new.data(), fun_diag_F.data(), Dimension::F, 1, Dimension::F);
            ARRAY_MATMUL(U_new.data(), Ubranch.data(), U.data(), Dimension::F, Dimension::F, Dimension::F);

            for (int i = 0; i < Dimension::F; ++i) c[i] = ((i == break_state) ? c[i] / norm_phase_b : 0.0e0);
            ARRAY_MATMUL_TRANS2(Ubranch.data(), c.data(), fun_diag_F.data(), Dimension::F, 1, Dimension::F);
            ARRAY_MATMUL(U.data(), Ubranch.data(), U.data(), Dimension::F, Dimension::F, Dimension::F);

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
        // PRINT_ARRAY(Snuc, Dimension::P, Dimension::P);
        // PRINT_ARRAY(c, Dimension::P, Dimension::F);
        // PRINT_ARRAY(Acoeff, Dimension::P, 1);
        // PRINT_ARRAY(rhored, Dimension::F, Dimension::F);
        // calc_density(rhored, Acoeff, Snuc, c, 1, 0, P_used, Dimension::P, Dimension::F);
        // PRINT_ARRAY(rhored, Dimension::F, Dimension::F);
    }

    // std::cout << "after P_used = " << P_used << "\n";
    return 0;

    // for (int iP = 0; iP < Dimension::P; ++iP) {
    //     span<kids_real> p    = this->p + iP * Dimension::N;
    //     span<kids_real> grad = this->grad + iP * Dimension::N;
    //     span<kids_real> dV   = this->dV + iP * Dimension::NFF;
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
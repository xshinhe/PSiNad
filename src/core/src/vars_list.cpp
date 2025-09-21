#include "psnd/vars_list.h"

// clang-format off

namespace PROJECT_NS {

std::vector<VARIABLE_BASE*> VARIABLE_BASE::_LIST;

namespace Dimension {
std::size_t M;      ///< Number of Monte Carlo calculations.
std::size_t P;      ///< Number of parallel trajectories (swarms of trajectories) in each Monte Carlo run.
std::size_t P_NOW;  ///< Active number of parallel trajectories (P_NOW <= P)
std::size_t N;      ///< Number of nuclear degrees of freedom. (=3*Natom for real molecular systems)
std::size_t F;      ///< Number of electronic degrees of freedom.
std::size_t N4;     ///< Full number of both nuclear & electronic DOFs (for monodromy matrix)
std::size_t Nb;     ///< Number of discretized modes for the bath
std::size_t nbath;  ///< Number of bathes coupled with systems
std::size_t L;      ///< Number of nonzero elements in interaction Q (intrinsic)
std::size_t MP;     ///< Product of M and P (M * P). @deprecated
std::size_t PP;     ///< Product of P and P (P * P). @deprecated
std::size_t PN;     ///< Product of P and N (P * N). @deprecated
std::size_t PNN;    ///< Product of P, N, and N (P * N * N). @deprecated
std::size_t PF;     ///< Product of P and F (P * F). @deprecated
std::size_t PFF;    ///< Product of P, F, and F (P * F * F). @deprecated
std::size_t PNF;    ///< Product of P, N, and F (P * N * F). @deprecated
std::size_t PNFF;   ///< Product of P, N, F, and F (P * N * F * F). @deprecated
std::size_t NF;     ///< Product of N and F (N * F). @deprecated
std::size_t NN;     ///< Product of N and N (N * N). @deprecated
std::size_t FF;     ///< Product of F and F (F * F). @deprecated
std::size_t NFF;    ///< Product of N, F, and F (N * F * F). @deprecated
std::size_t NNFF;   ///< Product of N, N, F, and F (N * N * F * F). @deprecated
std::size_t Fadd1;  ///< F plus 1 (F + 1). @deprecated
std::size_t sstep;  ///< Number of skipped steps for each sampling.
std::size_t nstep;  ///< Number of steps for a single simulation
std::size_t nsamp;  ///< Number of samplings for a single simulation
                                                    
Shape shape_1(1);                                  ///< Shape corresponding to a single element (1).
Shape shape_2(2);                                  ///< Shape corresponding to two elements (2).
Shape shape_X(1);                                  ///< Shape for an arbitrary number of elements.
Shape shape_Nxgrid(1);                             ///< Shape for Nxgrid. @deprecated
Shape shape_Nb(Shape::dynamic_t{&Nb});                       ///< Shape for the number of discretized modes.
Shape shape_nbathFF(Shape::dynamic_t{&nbath, &F, &F});       ///< Shape for the product of nbath, F, and F (nbath * F * F).
Shape shape_LNb(Shape::dynamic_t{&L, &Nb});                  ///< Shape for the product of L and Nb (L * Nb).
Shape shape_LnbathFF(Shape::dynamic_t{&L, &nbath, &F, &F});  ///< Shape for the product of L, nbath and F and F (L * nbath * F * F)
Shape shape_M(Shape::dynamic_t{&M});                         ///< Shape for the number of Monte Carlo calculations (M).
Shape shape_P(Shape::dynamic_t{&P});                         ///< Shape for the number of parallel trajectories (P).
Shape shape_N(Shape::dynamic_t{&N});                         ///< Shape for the number of nuclear degrees of freedom (N).
Shape shape_F(Shape::dynamic_t{&F});                         ///< Shape for the number of electronic degrees of freedom (F).
Shape shape_Fadd1(Shape::dynamic_t{&Fadd1});                 ///< Shape for F plus 1 (F + 1).
Shape shape_MP(Shape::dynamic_t{&M, &P});                    ///< Shape for the product of M and P (M * P).
Shape shape_PP(Shape::dynamic_t{&P, &P});                    ///< Shape for the product of P and P (P * P).
Shape shape_PN(Shape::dynamic_t{&P, &N});                    ///< Shape for the product of P and N (P * N).
Shape shape_PNN(Shape::dynamic_t{&P, &N, &N});               ///< Shape for the product of P, N, and N (P * N * N).
Shape shape_PF(Shape::dynamic_t{&P, &F});                    ///< Shape for the product of P and F (P * F).
Shape shape_PFF(Shape::dynamic_t{&P, &F, &F});               ///< Shape for the product of P, F, and F (P * F * F).
Shape shape_PNF(Shape::dynamic_t{&P, &N, &F});               ///< Shape for the product of P, N, and F (P * N * F).
Shape shape_PNFF(Shape::dynamic_t{&P, &N, &F, &F});          ///< Shape for the product of P, N, F, and F (P * N * F * F).
Shape shape_PN4N4(Shape::dynamic_t{&P, &N4, &N4});           ///< Shape for the product of (N + 2 * F)(N + 2 * F).
Shape shape_NF(Shape::dynamic_t{&N, &F});                    ///< Shape for the product of N and F (N * F).
Shape shape_NN(Shape::dynamic_t{&N, &N});                    ///< Shape for the product of N and N (N * N).
Shape shape_FF(Shape::dynamic_t{&F, &F});                    ///< Shape for the product of F and F (F * F).
Shape shape_NFF(Shape::dynamic_t{&N, &F, &F});               ///< Shape for the product of N, F, and F (N * F * F).
Shape shape_NNFF(Shape::dynamic_t{&N, &N, &F, &F});          ///< Shape for the product of N, N, F, and F (N * N * F * F).
Shape shape_PNNFF(Shape::dynamic_t{&P, &N, &N, &F, &F});     ///< Shape for the product of P, N, N, F, and F (P* N * N * F * F).

void static_build_shapes() {
    assert(M > 0 && P > 0 && N > 0 && F > 0);  //
    P_NOW = P;
    // auxiliary dimension
    FF    = F * F;
    NFF   = N * FF;
    NF    = N * F;
    N4    = N + 2 * F;
    NN    = N * N;
    PP    = P * P;
    PN    = P * N;
    PNN   = P * NN;
    PF    = P * F;
    PFF   = P * FF;
    PNF   = P * NF;
    PNFF  = P * NFF;
    MP    = M * P;
    Fadd1 = F + 1;
    shape_M.static_build();
    shape_P.static_build();
    shape_N.static_build();
    shape_F.static_build();
    shape_Fadd1.static_build();
    shape_MP.static_build();
    shape_PP.static_build();
    shape_PN.static_build();
    shape_PNN.static_build();
    shape_PF.static_build();
    shape_PFF.static_build();
    shape_PNFF.static_build();
    shape_NF.static_build();
    shape_NN.static_build();
    shape_FF.static_build();
    shape_NFF.static_build();
    shape_NNFF.static_build();
    shape_PNNFF.static_build();
    shape_PN4N4.static_build();
}
};  // namespace Dimension

namespace DATA {

#define NAME_WRAPPER(name, shape, doc_string) name(#name, shape, doc_string)
using namespace Dimension;
VARIABLE<psnd_real>    NAME_WRAPPER(init::Etot, &shape_P, "init total energy");
VARIABLE<psnd_real>    NAME_WRAPPER(init::p, &shape_PN, "init nuclear momentum");
VARIABLE<psnd_real>    NAME_WRAPPER(init::x, &shape_PN, "init nuclear coordinate");
VARIABLE<psnd_real>    NAME_WRAPPER(init::T, &shape_PFF, "init ADT matrix");
VARIABLE<psnd_complex> NAME_WRAPPER(init::c, &shape_PF, "init.c");
VARIABLE<psnd_complex> NAME_WRAPPER(init::cset, &shape_PFF, "init.cset");
VARIABLE<psnd_complex> NAME_WRAPPER(init::rho_ele, &shape_PFF, "init.rho_ele");
VARIABLE<psnd_complex> NAME_WRAPPER(init::rho_nuc, &shape_PFF, "init.rho_nuc");
VARIABLE<psnd_complex> NAME_WRAPPER(init::rho_dual, &shape_PFF, "init.rho_dual");

VARIABLE<psnd_real> NAME_WRAPPER(parameter::xi0, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::xi1, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::xi2, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::xi3, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::xiw, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::xir, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gamma0, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gamma1, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gamma2, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gamma3, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gammaw, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::gammar, &shape_1, "parameter");
VARIABLE<psnd_real> NAME_WRAPPER(parameter::Is, &shape_PFF, "parameter");

VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Acoeff, &shape_P, "configuration coefficients");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::E, &shape_PF, "adiabatic energies");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::Ekin, &shape_P, "kinematic energy");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::Epot, &shape_P, "potential energy");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::Etot, &shape_P, "total energy");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::GWP::L, &shape_P, "ES used for GWP");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::GWP::L1, &shape_P, "L1");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::GWP::L2, &shape_P, "L2");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::R, &shape_PP, "EV used for GWP");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::R1, &shape_PP, "R1");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::R2, &shape_PP, "R2");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::S1, &shape_PP, "overlap used for GWP");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::S1h, &shape_PP, "S1h");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::S2, &shape_PP, "S2");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::S2h, &shape_PP, "S2h");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::Sx, &shape_PP, "Sx");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::invS1h, &shape_PP, "inversed overlap used for GWP");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::GWP::invS2h, &shape_PP, "invS2h");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Hbasis, &shape_PP, "Hb");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Hcoeff, &shape_PP, "Hcoeff");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K0, &shape_PFF, "K0");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K1, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K1DA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K1DD, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K1QA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K1QD, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K2, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K2DA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K2DD, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K2QA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::K2QD, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::KSHA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::KTWA, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::KTWD, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::OpA, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::OpB, &shape_FF, "");
VARIABLE<psnd_int>     NAME_WRAPPER(integrator::P_used, &shape_1, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::S, &shape_PP, "S");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Sele, &shape_PP, "Sele");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Snuc, &shape_PP, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::U, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::UXdt, &shape_PP, "UX");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::UYdt, &shape_PP, "UY");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Ubranch, &shape_PFF, "Ub");  //? check size
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Udt, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Ucdt, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::Xcoeff, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::alpha, &shape_P, "");
// VARIABLE<psnd_real>    NAME_WRAPPER(integrator::alpha, &shape_N, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::c, &shape_PF, "electronic amplitude");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::cset, &shape_PFF, "electronic amplitude");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::c1, &shape_PN, "langevin coefficients 1");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::c2p, &shape_PN, "langevin coefficients 2");
VARIABLE<psnd_int>     NAME_WRAPPER(integrator::clone_account, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::dtAcoeff, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::dtSele, &shape_PP, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::dtlnSnuc, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::f, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::fadd, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::g, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::invS, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::m, &shape_PN, "mass of trajectories");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::minv, &shape_PN, "inversed mass of trajectories");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::nhc::G, &shape_X, "G");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::nhc::Q, &shape_X, "Q");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::nhc::p, &shape_X, "p");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::nhc::x, &shape_X, "x");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::norm, &shape_P, "");
VARIABLE<psnd_int>     NAME_WRAPPER(integrator::occ_nuc, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::ve, &shape_PN, "velocity of trajectories");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::p, &shape_PN, "momentum of trajectories");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::p_sign, &shape_2, "");  // check size? 2P?
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::param::c1, &shape_X, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::param::c2p, &shape_X, "");
VARIABLE<psnd_bint>    NAME_WRAPPER(integrator::pf_cross, &shape_PF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rho_dual, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rho_ele, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rho_nuc, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rhored, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rhored2, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::rhored3, &shape_FF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::trKTWA, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::trKTWD, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::sqcw, &shape_PF, "");
// VARIABLE<psnd_real>    NAME_WRAPPER(integrator::sqcw, &shape_FF, "");
// VARIABLE<psnd_real>    NAME_WRAPPER(integrator::sqcw, &shape_F, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::sqcw0, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::sqcwh, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::monodromy::mono, &shape_PN4N4, "mono");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::monodromy::monodt, &shape_PN4N4, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp1, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp2, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp3, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp4, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp5, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::monodromy::MFFtmp6, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::forceeval::mask, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::forceeval::dmask, &shape_PNFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::I_PP, &shape_PP, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::MatC_PP, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::MatR_PP, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::TtTold, &shape_FF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::direction, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::fproj, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::ftmp, &shape_N, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::fun_diag_F, &shape_F, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::fun_diag_P, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::invexpidiagdt, &shape_F, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::ve, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::tmp::vedE, &shape_FF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::tmp::wrho, &shape_FF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::relwgt, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::gf_x, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::gf_p, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::gf_c, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::gf_all, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::avgx, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::varx, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::avgp, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::varp, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::avgxf, &shape_NF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::varxf, &shape_NF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::xintercept, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::xinterceptf, &shape_NFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::xslope, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::term_1, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::term_2, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::fadiat, &shape_PNF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::COUP::pb, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::veF, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_AA, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_AD, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_CC, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_CP, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_DD, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::w_PP, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::ww_A, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::ww_D, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::wz_A, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(integrator::wz_D, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::x, &shape_PN, "coordinate of trajectories");
VARIABLE<psnd_real>    NAME_WRAPPER(integrator::xgrid, &shape_Nxgrid, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::scheme_id, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::solver_id, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::calc_id, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::step_id, &shape_1, "");
VARIABLE<psnd_bint>    NAME_WRAPPER(control::at_condition, &shape_1, "");
// VARIABLE<psnd_bint>    NAME_WRAPPER(control::at_samplingstep_finally, &shape_1, "");
// VARIABLE<psnd_bint>    NAME_WRAPPER(control::at_samplingstep_initially, &shape_1, "");
VARIABLE<psnd_real>    NAME_WRAPPER(control::dt, &shape_1, "");
VARIABLE<psnd_real>    NAME_WRAPPER(control::pertimeunit, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::dtsize, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::last_tried_dtsize, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::fail_type, &shape_1, "");
VARIABLE<psnd_bint>    NAME_WRAPPER(control::frez, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::isamp, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::istep, &shape_1, "");
VARIABLE<psnd_bint>    NAME_WRAPPER(control::last_attempt, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::msize, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::nsamp, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::nstep, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::sstep, &shape_1, "");
VARIABLE<psnd_bint>    NAME_WRAPPER(control::succ, &shape_1, "");
VARIABLE<psnd_real>    NAME_WRAPPER(control::t, &shape_1, "");
VARIABLE<psnd_int>     NAME_WRAPPER(control::tsize, &shape_1, "");
VARIABLE<psnd_real>    NAME_WRAPPER(last::Etot, &shape_P, "total energy in last step");
VARIABLE<psnd_complex> NAME_WRAPPER(last::c, &shape_PF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(last::dV, &shape_PNFF, "dV");
VARIABLE<psnd_real>    NAME_WRAPPER(last::g, &shape_P, "");
VARIABLE<psnd_real>    NAME_WRAPPER(last::grad, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(last::p, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(last::x, &shape_PN, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::Hsys, &shape_FF, "system part of Hamiltonian");
VARIABLE<psnd_real>    NAME_WRAPPER(model::Kmat, &shape_NN, "nuclear oscillation matrix");
VARIABLE<psnd_real>    NAME_WRAPPER(model::Qmat, &shape_NFF, "coupling matrix");
VARIABLE<psnd_real>    NAME_WRAPPER(model::Tmod, &shape_NN, "normalmode transformation");
VARIABLE<psnd_real>    NAME_WRAPPER(model::V, &shape_PFF, "diabatic potential of energy matrix");
VARIABLE<psnd_int>     NAME_WRAPPER(model::atoms, &shape_N, "atom indices");          // only N//3 size
VARIABLE<psnd_int>     NAME_WRAPPER(model::layer_type, &shape_N, "atom layer type");  // only N//3 size
VARIABLE<psnd_real>    NAME_WRAPPER(model::bath::coeffs, &shape_Nb, "discretized coefficients for a bath");
VARIABLE<psnd_real>    NAME_WRAPPER(model::bath::omegas, &shape_Nb, "discretized frequencies for a bath");
VARIABLE<psnd_real>    NAME_WRAPPER(model::coupling::CL, &shape_LNb, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::coupling::Q, &shape_nbathFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::coupling::QL, &shape_LnbathFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::dV, &shape_PNFF, "dV");
VARIABLE<psnd_real>    NAME_WRAPPER(model::ddV, &shape_PNNFF, "ddV");
VARIABLE<psnd_real>    NAME_WRAPPER(model::f_p, &shape_F, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::f_r, &shape_F, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::f_rp, &shape_F, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::grad, &shape_PN, "grad");
VARIABLE<psnd_real>    NAME_WRAPPER(model::hess, &shape_PNN, "hess");
VARIABLE<psnd_real>    NAME_WRAPPER(model::kcoeff, &shape_X, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::lcoeff, &shape_X, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::mass, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::p0, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::p_sigma, &shape_N, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::Vc, &shape_PFF, "Vc");
VARIABLE<psnd_complex> NAME_WRAPPER(model::dVc, &shape_PNFF, "dVc");
VARIABLE<psnd_complex> NAME_WRAPPER(model::ddVc, &shape_PNNFF, "ddVc");
VARIABLE<psnd_real>    NAME_WRAPPER(model::MB::Jpmat, &shape_PP, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::MB::Jzmat, &shape_PP, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::MB::SXred, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::MB::SYred, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::MB::SZred, &shape_P, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::MB::H1, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::MB::H2, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::eig, &shape_PF, "eig");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::E, &shape_PFF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::rep::H, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::lam, &shape_PF, "");
VARIABLE<psnd_complex> NAME_WRAPPER(model::rep::R, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::T, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::Told, &shape_PFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::dE, &shape_PNFF, "dE");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::ddE, &shape_PNNFF, "ddE");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::nac, &shape_PNFF, "nac");
VARIABLE<psnd_real>    NAME_WRAPPER(model::rep::nac_prev, &shape_PNFF, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::vpes, &shape_P, "scalar potential of energy");
VARIABLE<psnd_real>    NAME_WRAPPER(model::w, &shape_N, "");
VARIABLE<psnd_real>    NAME_WRAPPER(model::x0, &shape_N, "initial cencter of coordinate");
VARIABLE<psnd_real>    NAME_WRAPPER(model::x_sigma, &shape_N, "initial width of coordinate");
VARIABLE<psnd_int>     NAME_WRAPPER(random::seed, &shape_1, "");

// custom variables


};  // namespace DATA

};  // namespace PROJECT_NS

// clang-format on

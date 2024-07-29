/**@file        vars_list.h
 * @brief       declaration of variables used in the program.
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> Initial version.
 * </table>
 **********************************************************************************
 */


#ifndef VARS_LIST_H
#define VARS_LIST_H

#include "kids/Shape.h"
#include "kids/Types.h"
#include "kids/Variable.h"

namespace PROJECT_NS {

/**
 * Dimension namespace locates read parameters for system size
 */
namespace Dimension {
extern std::size_t M;  ///< Number of Monte Carlo calculations.
extern std::size_t P;  ///< Number of parallel trajectories (swarms of trajectories) in each Monte Carlo run.
extern std::size_t P_NOW;
extern std::size_t N;      ///< Number of nuclear degrees of freedom.
extern std::size_t F;      ///< Number of electronic degrees of freedom.
extern std::size_t Nb;     ///< Number of discretized modes
extern std::size_t nbath;  ///< Number of bathes
extern std::size_t L;      ///< Number of nonzero elements in interaction Q

extern std::size_t MP;    ///< Product of M and P (M * P).
extern std::size_t PP;    ///< Product of P and P (P * P).
extern std::size_t PN;    ///< Product of P and N (P * N).
extern std::size_t PNN;   ///< Product of P, N, and N (P * N * N).
extern std::size_t PF;    ///< Product of P and F (P * F).
extern std::size_t PFF;   ///< Product of P, F, and F (P * F * F).
extern std::size_t PNFF;  ///< Product of P, N, F, and F (P * N * F * F).
extern std::size_t NF;    ///< Product of N and F (N * F).
extern std::size_t NN;    ///< Product of N and N (N * N).
extern std::size_t FF;    ///< Product of F and F (F * F).
extern std::size_t NFF;   ///< Product of N, F, and F (N * F * F).
extern std::size_t NNFF;  ///< Product of N, N, F, and F (N * N * F * F).
extern std::size_t N4;    ///<  (N + 2* F).

extern std::size_t Fadd1;  ///< F plus 1 (F + 1).

extern std::size_t sstep;
extern std::size_t nstep;
extern std::size_t nsamp;

extern Shape shape_1;      ///< Shape corresponding to a single element (1).
extern Shape shape_2;      ///< Shape corresponding to two elements (2).
extern Shape shape_X;      ///< Shape for an arbitrary number of elements.
extern Shape shape_M;      ///< Shape for the number of Monte Carlo calculations (M).
extern Shape shape_P;      ///< Shape for the number of parallel trajectories (P).
extern Shape shape_N;      ///< Shape for the number of nuclear degrees of freedom (N).
extern Shape shape_F;      ///< Shape for the number of electronic degrees of freedom (F).
extern Shape shape_Fadd1;  ///< Shape for F plus 1 (F + 1).
extern Shape shape_MP;     ///< Shape for the product of M and P (M * P).
extern Shape shape_PP;     ///< Shape for the product of P and P (P * P).
extern Shape shape_PN;     ///< Shape for the product of P and N (P * N).
extern Shape shape_PNN;    ///< Shape for the product of P, N, and N (P * N * N).
extern Shape shape_PF;     ///< Shape for the product of P and F (P * F).
extern Shape shape_PFF;    ///< Shape for the product of P, F, and F (P * F * F).
extern Shape shape_PNFF;   ///< Shape for the product of P, N, F, and F (P * N * F * F).
extern Shape shape_NF;     ///< Shape for the product of N and F (N * F).
extern Shape shape_NN;     ///< Shape for the product of N and N (N * N).
extern Shape shape_FF;     ///< Shape for the product of F and F (F * F).
extern Shape shape_NFF;    ///< Shape for the product of N, F, and F (N * F * F).
extern Shape shape_NNFF;   ///< Shape for the product of N, N, F, and F (N * N * F * F).
extern Shape shape_PN4N4;  ///< Shape for the product of (N + 2 * F)(N + 2 * F).

extern Shape shape_Nb;        ///< Shape for the number of discretized modes
extern Shape shape_nbathFF;   ///< Shape for the product of nbath, F, and F (nbath * F * F)
extern Shape shape_LNb;       ///< Shape for the product of L and Nb (L * Nb)
extern Shape shape_LnbathFF;  ///< Shape for the product of L, nbath and F and F (L * nbath * F * F)

extern void static_build_shapes();
};  // namespace Dimension

namespace DATA {

namespace init {  // @to be removed
extern VARIABLE<kids_real>    Etot;
extern VARIABLE<kids_real>    p;
extern VARIABLE<kids_real>    x;
extern VARIABLE<kids_real>    T;
extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_complex> cset;
extern VARIABLE<kids_complex> rho_ele;
extern VARIABLE<kids_complex> rho_nuc;
extern VARIABLE<kids_complex> rho_dual;
};  // namespace init

namespace parameter {
extern VARIABLE<kids_real> xi0;
extern VARIABLE<kids_real> xi1;
extern VARIABLE<kids_real> xi2;
extern VARIABLE<kids_real> xi3;
extern VARIABLE<kids_real> xiw;
extern VARIABLE<kids_real> xir;
extern VARIABLE<kids_real> gamma0;
extern VARIABLE<kids_real> gamma1;
extern VARIABLE<kids_real> gamma2;
extern VARIABLE<kids_real> gamma3;
extern VARIABLE<kids_real> gammaw;
extern VARIABLE<kids_real> gammar;
extern VARIABLE<kids_real> Is;
};  // namespace parameter

namespace integrator {
extern VARIABLE<kids_complex> Acoeff;
extern VARIABLE<kids_real>    E;
extern VARIABLE<kids_real>    Ekin;
extern VARIABLE<kids_real>    Epot;
extern VARIABLE<kids_real>    Etot;

namespace GWP {
extern VARIABLE<kids_real>    L;
extern VARIABLE<kids_real>    L1;
extern VARIABLE<kids_real>    L2;
extern VARIABLE<kids_complex> R;
extern VARIABLE<kids_complex> R1;
extern VARIABLE<kids_complex> R2;
extern VARIABLE<kids_complex> S1;
extern VARIABLE<kids_complex> S1h;
extern VARIABLE<kids_complex> S2;
extern VARIABLE<kids_complex> S2h;
extern VARIABLE<kids_complex> Sx;
extern VARIABLE<kids_complex> invS1h;
extern VARIABLE<kids_complex> invS2h;
};  // namespace GWP

extern VARIABLE<kids_complex> Hbasis;
extern VARIABLE<kids_complex> Hcoeff;
extern VARIABLE<kids_complex> K0;
extern VARIABLE<kids_complex> K1;
extern VARIABLE<kids_complex> K1DA;
extern VARIABLE<kids_complex> K1DD;
extern VARIABLE<kids_complex> K1QA;
extern VARIABLE<kids_complex> K1QD;
extern VARIABLE<kids_complex> K2;
extern VARIABLE<kids_complex> K2DA;
extern VARIABLE<kids_complex> K2DD;
extern VARIABLE<kids_complex> K2QA;
extern VARIABLE<kids_complex> K2QD;
extern VARIABLE<kids_complex> KSHA;
extern VARIABLE<kids_complex> KTWA;
extern VARIABLE<kids_complex> KTWD;
extern VARIABLE<kids_complex> OpA;
extern VARIABLE<kids_complex> OpB;
extern VARIABLE<kids_real>    P_used;
extern VARIABLE<kids_complex> S;
extern VARIABLE<kids_complex> Sele;
extern VARIABLE<kids_complex> Snuc;
extern VARIABLE<kids_complex> U;
extern VARIABLE<kids_complex> UXdt;
extern VARIABLE<kids_complex> UYdt;
extern VARIABLE<kids_complex> Ubranch;
extern VARIABLE<kids_complex> Udt;
extern VARIABLE<kids_complex> Xcoeff;
extern VARIABLE<kids_real>    alpha;
extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_complex> cset;
extern VARIABLE<kids_real>    c1;
extern VARIABLE<kids_real>    c2p;
extern VARIABLE<kids_int>     clone_account;
extern VARIABLE<kids_complex> dtAcoeff;
extern VARIABLE<kids_complex> dtSele;
extern VARIABLE<kids_complex> dtlnSnuc;
extern VARIABLE<kids_real>    f;
extern VARIABLE<kids_real>    fadd;
extern VARIABLE<kids_real>    g;
extern VARIABLE<kids_complex> invS;
extern VARIABLE<kids_real>    m;
extern VARIABLE<kids_real>    minv;

namespace nhc {
extern VARIABLE<kids_real> G;
extern VARIABLE<kids_real> Q;
extern VARIABLE<kids_real> p;
extern VARIABLE<kids_real> x;
};  // namespace nhc

extern VARIABLE<kids_real>    norm;
extern VARIABLE<kids_int>     occ_nuc;
extern VARIABLE<kids_real>    ve;
extern VARIABLE<kids_real>    p;
extern VARIABLE<kids_complex> p_sign;

namespace param {
extern VARIABLE<kids_real> c1;
extern VARIABLE<kids_real> c2p;
};  // namespace param

extern VARIABLE<kids_bint>    pf_cross;
extern VARIABLE<kids_complex> rho_dual;
extern VARIABLE<kids_complex> rho_ele;
extern VARIABLE<kids_complex> rho_nuc;
extern VARIABLE<kids_complex> rhored;
extern VARIABLE<kids_complex> rhored2;
extern VARIABLE<kids_complex> rhored3;
extern VARIABLE<kids_real>    trKTWA;
extern VARIABLE<kids_real>    trKTWD;
extern VARIABLE<kids_real>    sqcw;
extern VARIABLE<kids_real>    sqcw0;
extern VARIABLE<kids_real>    sqcwh;

namespace monodromy {
extern VARIABLE<kids_real>    mono;
extern VARIABLE<kids_real>    monodt;
extern VARIABLE<kids_complex> MFFtmp1;
extern VARIABLE<kids_complex> MFFtmp2;
extern VARIABLE<kids_complex> MFFtmp3;
extern VARIABLE<kids_complex> MFFtmp4;
extern VARIABLE<kids_complex> MFFtmp5;
extern VARIABLE<kids_complex> MFFtmp6;
};  // namespace monodromy

namespace forceeval {
extern VARIABLE<kids_complex> mask;
extern VARIABLE<kids_complex> dmask;
};  // namespace forceeval

namespace tmp {
extern VARIABLE<kids_complex> I_PP;
extern VARIABLE<kids_complex> MatC_PP;
extern VARIABLE<kids_real>    MatR_PP;
extern VARIABLE<kids_real>    TtTold;
extern VARIABLE<kids_real>    direction;
extern VARIABLE<kids_real>    fproj;
extern VARIABLE<kids_real>    ftmp;
extern VARIABLE<kids_complex> fun_diag_F;
extern VARIABLE<kids_complex> fun_diag_P;
extern VARIABLE<kids_complex> invexpidiagdt;
extern VARIABLE<kids_real>    ve;
extern VARIABLE<kids_real>    vedE;
extern VARIABLE<kids_complex> wrho;
};  // namespace tmp

extern VARIABLE<kids_real>    veF;
extern VARIABLE<kids_complex> w;
extern VARIABLE<kids_complex> w_AA;
extern VARIABLE<kids_complex> w_AD;
extern VARIABLE<kids_complex> w_CC;
extern VARIABLE<kids_complex> w_CP;
extern VARIABLE<kids_complex> w_DD;
extern VARIABLE<kids_complex> w_PP;
extern VARIABLE<kids_complex> ww_A;
extern VARIABLE<kids_complex> ww_D;
extern VARIABLE<kids_complex> wz_A;
extern VARIABLE<kids_complex> wz_D;
extern VARIABLE<kids_real>    x;
extern VARIABLE<kids_real>    xgrid;
};  // namespace integrator


namespace flowcontrol {
extern VARIABLE<kids_int>  scheme_id;
extern VARIABLE<kids_int>  solver_id;
extern VARIABLE<kids_int>  calc_id;
extern VARIABLE<kids_int>  step_id;
extern VARIABLE<kids_bint> at_condition;
// extern VARIABLE<kids_bint> at_samplingstep_finally;
// extern VARIABLE<kids_bint> at_samplingstep_initially;
extern VARIABLE<kids_real> dt;
extern VARIABLE<kids_real> pertimeunit;
extern VARIABLE<kids_int>  dtsize;
extern VARIABLE<kids_int>  fail_type;
extern VARIABLE<kids_bint> frez;
extern VARIABLE<kids_int>  isamp;
extern VARIABLE<kids_int>  istep;
extern VARIABLE<kids_bint> last_attempt;
extern VARIABLE<kids_int>  msize;
extern VARIABLE<kids_int>  nsamp;
extern VARIABLE<kids_int>  nstep;
extern VARIABLE<kids_int>  sstep;
extern VARIABLE<kids_bint> succ;
extern VARIABLE<kids_real> t;
extern VARIABLE<kids_int>  tsize;
};  // namespace flowcontrol


namespace last {
extern VARIABLE<kids_real>    Etot;
extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_real>    dV;
extern VARIABLE<kids_real>    g;
extern VARIABLE<kids_real>    grad;
extern VARIABLE<kids_real>    p;
extern VARIABLE<kids_real>    x;
};  // namespace last


namespace model {
extern VARIABLE<kids_real> Hsys;
extern VARIABLE<kids_real> Kmat;
extern VARIABLE<kids_real> Qmat;
extern VARIABLE<kids_real> Tmod;
extern VARIABLE<kids_real> V;
extern VARIABLE<kids_int>  atoms;

namespace bath {
extern VARIABLE<kids_real> coeffs;
extern VARIABLE<kids_real> omegas;
};  // namespace bath


namespace coupling {
extern VARIABLE<kids_real> CL;
extern VARIABLE<kids_real> Q;
extern VARIABLE<kids_real> QL;
};  // namespace coupling

extern VARIABLE<kids_real> dV;
extern VARIABLE<kids_real> ddV;
extern VARIABLE<kids_real> f_p;
extern VARIABLE<kids_real> f_r;
extern VARIABLE<kids_real> f_rp;
extern VARIABLE<kids_real> grad;
extern VARIABLE<kids_real> hess;
extern VARIABLE<kids_real> kcoeff;
extern VARIABLE<kids_real> lcoeff;
extern VARIABLE<kids_real> mass;
extern VARIABLE<kids_real> p0;
extern VARIABLE<kids_real> p_sigma;

namespace MB {
extern VARIABLE<kids_real>    Jpmat;
extern VARIABLE<kids_real>    Jzmat;
extern VARIABLE<kids_complex> SXred;
extern VARIABLE<kids_complex> SYred;
extern VARIABLE<kids_complex> SZred;
extern VARIABLE<kids_complex> H1;
extern VARIABLE<kids_complex> H2;
};  // namespace MB

namespace rep {
extern VARIABLE<kids_real>    eig;
extern VARIABLE<kids_real>    E;
extern VARIABLE<kids_complex> H;
extern VARIABLE<kids_real>    lam;
extern VARIABLE<kids_complex> R;
extern VARIABLE<kids_real>    T;
extern VARIABLE<kids_real>    Told;
extern VARIABLE<kids_real>    dE;
extern VARIABLE<kids_real>    ddE;
extern VARIABLE<kids_real>    nac;
extern VARIABLE<kids_real>    nac_prev;
};  // namespace rep

extern VARIABLE<kids_real> vpes;
extern VARIABLE<kids_real> w;
extern VARIABLE<kids_real> x0;
extern VARIABLE<kids_real> x_sigma;
};  // namespace model


namespace random {
extern VARIABLE<kids_int> seed;
};  // namespace random

};  // namespace DATA

};  // namespace PROJECT_NS

#endif  // VARS_LIST_H

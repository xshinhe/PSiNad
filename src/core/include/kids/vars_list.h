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

namespace Dimension {
extern std::size_t M;  ///< Number of Monte Carlo calculations.
};  // namespace Dimension
namespace Dimension {
extern std::size_t P;  ///< Number of parallel trajectories (swarms of trajectories) in each Monte Carlo run.
};  // namespace Dimension
namespace Dimension {
extern std::size_t P_NOW;  ///< Active number of parallel trajectories (P_NOW <= P)
};  // namespace Dimension
namespace Dimension {
extern std::size_t N;  ///< Number of nuclear degrees of freedom. (=3*Natom for real molecular systems)
};  // namespace Dimension
namespace Dimension {
extern std::size_t F;  ///< Number of electronic degrees of freedom.
};  // namespace Dimension
namespace Dimension {
extern std::size_t N4;  ///< Full number of both nuclear & electronic DOFs (for monodromy matrix)
};  // namespace Dimension
namespace Dimension {
extern std::size_t Nb;  ///< Number of discretized modes for the bath
};  // namespace Dimension
namespace Dimension {
extern std::size_t nbath;  ///< Number of bathes coupled with systems
};  // namespace Dimension
namespace Dimension {
extern std::size_t L;  ///< Number of nonzero elements in interaction Q (intrinsic)
};  // namespace Dimension
namespace Dimension {
extern std::size_t MP;  ///< Product of M and P (M * P). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PP;  ///< Product of P and P (P * P). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PN;  ///< Product of P and N (P * N). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PNN;  ///< Product of P, N, and N (P * N * N). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PF;  ///< Product of P and F (P * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PFF;  ///< Product of P, F, and F (P * F * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PNF;  ///< Product of P, N, and F (P * N * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t PNFF;  ///< Product of P, N, F, and F (P * N * F * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t NF;  ///< Product of N and F (N * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t NN;  ///< Product of N and N (N * N). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t FF;  ///< Product of F and F (F * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t NFF;  ///< Product of N, F, and F (N * F * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t NNFF;  ///< Product of N, N, F, and F (N * N * F * F). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t Fadd1;  ///< F plus 1 (F + 1). @deprecated
};  // namespace Dimension
namespace Dimension {
extern std::size_t sstep;  ///< Number of skipped steps for each sampling.
};  // namespace Dimension
namespace Dimension {
extern std::size_t nstep;  ///< Number of steps for a single simulation
};  // namespace Dimension
namespace Dimension {
extern std::size_t nsamp;  ///< Number of samplings for a single simulation
};  // namespace Dimension
namespace Dimension {
extern Shape shape_1;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_2;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_X;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_Nxgrid;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_Nb;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_nbathFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_LNb;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_LnbathFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_M;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_P;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_N;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_F;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_Fadd1;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_MP;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PP;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PN;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PNN;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PNF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PNFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PN4N4;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_NF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_NN;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_FF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_NFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_NNFF;
};  // namespace Dimension
namespace Dimension {
extern Shape shape_PNNFF;
};  // namespace Dimension
namespace DATA {
namespace init {
extern VARIABLE<kids_real> Etot;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_real> p;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_real> x;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_real> T;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_complex> c;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_complex> cset;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_complex> rho_ele;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_complex> rho_nuc;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<kids_complex> rho_dual;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xi0;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xi1;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xi2;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xi3;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xiw;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> xir;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gamma0;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gamma1;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gamma2;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gamma3;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gammaw;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> gammar;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<kids_real> Is;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Acoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> E;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> Ekin;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> Epot;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> Etot;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_real> L;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_real> L1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_real> L2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> R;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> R1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> R2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> S1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> S1h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> S2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> S2h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> Sx;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> invS1h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<kids_complex> invS2h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Hbasis;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Hcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K0;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K1;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K1DA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K1DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K1QA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K1QD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K2;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K2DA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K2DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K2QA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> K2QD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> KSHA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> KTWA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> KTWD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> OpA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> OpB;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_int> P_used;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> S;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Sele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Snuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> U;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> UXdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> UYdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Ubranch;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Udt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Ucdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> Xcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> alpha;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> c;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> cset;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> c1;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> c2p;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_int> clone_account;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> dtAcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> dtSele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> dtlnSnuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> f;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> fadd;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> g;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> invS;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> m;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> minv;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<kids_real> G;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<kids_real> Q;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<kids_real> p;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<kids_real> x;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> norm;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_int> occ_nuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> ve;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> p;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> p_sign;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace param {
extern VARIABLE<kids_real> c1;
};  // namespace param
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace param {
extern VARIABLE<kids_real> c2p;
};  // namespace param
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_bint> pf_cross;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rho_dual;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rho_ele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rho_nuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rhored;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rhored2;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> rhored3;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> trKTWA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> trKTWD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> sqcw;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> sqcw0;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> sqcwh;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_real> mono;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_real> monodt;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp1;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp2;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp3;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp4;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp5;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<kids_complex> MFFtmp6;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace forceeval {
extern VARIABLE<kids_complex> mask;
};  // namespace forceeval
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace forceeval {
extern VARIABLE<kids_complex> dmask;
};  // namespace forceeval
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> I_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> MatC_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> MatR_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> TtTold;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> direction;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> fproj;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> ftmp;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> fun_diag_F;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> fun_diag_P;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> invexpidiagdt;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> ve;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_real> vedE;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<kids_complex> wrho;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> relwgt;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> gf_x;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> gf_p;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> gf_c;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> gf_all;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> avgx;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> varx;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> avgp;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> varp;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> avgxf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> varxf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> xintercept;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> xinterceptf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> xslope;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> term_1;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> term_2;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> fadiat;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<kids_real> pb;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> veF;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_AA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_AD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_CC;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_CP;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> w_PP;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> ww_A;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> ww_D;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> wz_A;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_complex> wz_D;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> x;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<kids_real> xgrid;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> scheme_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> solver_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> calc_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> step_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_bint> at_condition;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_real> dt;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_real> pertimeunit;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> dtsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> last_tried_dtsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> fail_type;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_bint> frez;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> isamp;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> istep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_bint> last_attempt;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> msize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> nsamp;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> nstep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> sstep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_bint> succ;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_real> t;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<kids_int> tsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> Etot;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_complex> c;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> dV;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> g;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> grad;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> p;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<kids_real> x;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> Hsys;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> Kmat;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> Qmat;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> Tmod;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> V;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_int> atoms;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_int> layer_type;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace bath {
extern VARIABLE<kids_real> coeffs;
};  // namespace bath
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace bath {
extern VARIABLE<kids_real> omegas;
};  // namespace bath
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<kids_real> CL;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<kids_real> Q;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<kids_real> QL;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> dV;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> ddV;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> f_p;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> f_r;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> f_rp;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> grad;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> hess;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> kcoeff;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> lcoeff;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> mass;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> p0;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> p_sigma;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_complex> Vc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_complex> dVc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_complex> ddVc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_real> Jpmat;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_real> Jzmat;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_complex> SXred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_complex> SYred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_complex> SZred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_complex> H1;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<kids_complex> H2;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> eig;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> E;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_complex> H;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> lam;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_complex> R;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> T;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> Told;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> dE;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> ddE;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> nac;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<kids_real> nac_prev;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> vpes;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> w;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> x0;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<kids_real> x_sigma;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace random {
extern VARIABLE<kids_int> seed;
};  // namespace random
};  // namespace DATA
namespace Dimension {
extern void static_build_shapes();
};  // namespace Dimension
};  // namespace PROJECT_NS
#endif  // VARS_LIST_H

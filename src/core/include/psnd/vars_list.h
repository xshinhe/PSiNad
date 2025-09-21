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
 *  This software is a product of academic research conducted by Professor Liu's
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

#ifndef PSND_VARS_LIST_H
#define PSND_VARS_LIST_H

#include "psnd/Shape.h"
#include "psnd/Types.h"
#include "psnd/Variable.h"

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
extern VARIABLE<psnd_real> Etot;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_real> p;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_real> x;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_real> T;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_complex> c;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_complex> cset;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_complex> rho_ele;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_complex> rho_nuc;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace init {
extern VARIABLE<psnd_complex> rho_dual;
};  // namespace init
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xi0;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xi1;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xi2;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xi3;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xiw;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> xir;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gamma0;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gamma1;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gamma2;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gamma3;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gammaw;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> gammar;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace parameter {
extern VARIABLE<psnd_real> Is;
};  // namespace parameter
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Acoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> E;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> Ekin;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> Epot;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> Etot;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_real> L;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_real> L1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_real> L2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> R;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> R1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> R2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> S1;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> S1h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> S2;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> S2h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> Sx;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> invS1h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace GWP {
extern VARIABLE<psnd_complex> invS2h;
};  // namespace GWP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Hbasis;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Hcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K0;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K1;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K1DA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K1DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K1QA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K1QD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K2;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K2DA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K2DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K2QA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> K2QD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> KSHA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> KTWA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> KTWD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> OpA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> OpB;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_int> P_used;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> S;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Sele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Snuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> U;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> UXdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> UYdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Ubranch;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Udt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Ucdt;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> Xcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> alpha;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> c;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> cset;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> c1;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> c2p;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_int> clone_account;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> dtAcoeff;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> dtSele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> dtlnSnuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> f;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> fadd;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> g;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> invS;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> m;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> minv;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<psnd_real> G;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<psnd_real> Q;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<psnd_real> p;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace nhc {
extern VARIABLE<psnd_real> x;
};  // namespace nhc
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> norm;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_int> occ_nuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> ve;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> p;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> p_sign;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace param {
extern VARIABLE<psnd_real> c1;
};  // namespace param
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace param {
extern VARIABLE<psnd_real> c2p;
};  // namespace param
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_bint> pf_cross;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rho_dual;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rho_ele;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rho_nuc;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rhored;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rhored2;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> rhored3;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> trKTWA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> trKTWD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> sqcw;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> sqcw0;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> sqcwh;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_real> mono;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_real> monodt;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp1;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp2;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp3;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp4;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp5;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace monodromy {
extern VARIABLE<psnd_complex> MFFtmp6;
};  // namespace monodromy
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace forceeval {
extern VARIABLE<psnd_complex> mask;
};  // namespace forceeval
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace forceeval {
extern VARIABLE<psnd_complex> dmask;
};  // namespace forceeval
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> I_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> MatC_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> MatR_PP;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> TtTold;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> direction;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> fproj;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> ftmp;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> fun_diag_F;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> fun_diag_P;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> invexpidiagdt;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> ve;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_real> vedE;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace tmp {
extern VARIABLE<psnd_complex> wrho;
};  // namespace tmp
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> relwgt;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> gf_x;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> gf_p;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> gf_c;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> gf_all;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> avgx;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> varx;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> avgp;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> varp;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> avgxf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> varxf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> xintercept;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> xinterceptf;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> xslope;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> term_1;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> term_2;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> fadiat;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
namespace COUP {
extern VARIABLE<psnd_real> pb;
};  // namespace COUP
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> veF;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_AA;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_AD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_CC;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_CP;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_DD;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> w_PP;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> ww_A;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> ww_D;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> wz_A;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_complex> wz_D;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> x;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace integrator {
extern VARIABLE<psnd_real> xgrid;
};  // namespace integrator
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> scheme_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> solver_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> calc_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> step_id;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_bint> at_condition;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_real> dt;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_real> pertimeunit;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> dtsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> last_tried_dtsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> fail_type;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_bint> frez;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> isamp;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> istep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_bint> last_attempt;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> msize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> nsamp;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> nstep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> sstep;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_bint> succ;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_real> t;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace control {
extern VARIABLE<psnd_int> tsize;
};  // namespace control
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> Etot;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_complex> c;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> dV;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> g;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> grad;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> p;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace last {
extern VARIABLE<psnd_real> x;
};  // namespace last
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> Hsys;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> Kmat;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> Qmat;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> Tmod;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> V;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_int> atoms;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_int> layer_type;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace bath {
extern VARIABLE<psnd_real> coeffs;
};  // namespace bath
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace bath {
extern VARIABLE<psnd_real> omegas;
};  // namespace bath
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<psnd_real> CL;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<psnd_real> Q;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace coupling {
extern VARIABLE<psnd_real> QL;
};  // namespace coupling
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> dV;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> ddV;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> f_p;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> f_r;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> f_rp;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> grad;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> hess;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> kcoeff;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> lcoeff;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> mass;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> p0;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> p_sigma;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_complex> Vc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_complex> dVc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_complex> ddVc;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_real> Jpmat;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_real> Jzmat;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_complex> SXred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_complex> SYred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_complex> SZred;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_complex> H1;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace MB {
extern VARIABLE<psnd_complex> H2;
};  // namespace MB
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> eig;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> E;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_complex> H;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> lam;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_complex> R;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> T;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> Told;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> dE;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> ddE;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> nac;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
namespace rep {
extern VARIABLE<psnd_real> nac_prev;
};  // namespace rep
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> vpes;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> w;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> x0;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace model {
extern VARIABLE<psnd_real> x_sigma;
};  // namespace model
};  // namespace DATA
namespace DATA {
namespace random {
extern VARIABLE<psnd_int> seed;
};  // namespace random
};  // namespace DATA
namespace Dimension {
extern void static_build_shapes();
};  // namespace Dimension
};  // namespace PROJECT_NS
#endif  // PSND_VARS_LIST_H

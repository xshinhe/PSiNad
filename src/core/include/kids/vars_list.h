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
extern std::size_t P;  ///< Number of parallel trajectories in each run
extern std::size_t N;  ///< Number of nuclear degrees of freedom.
extern std::size_t F;  ///< Number of electronic degrees of freedom.
extern std::size_t MP;
extern std::size_t PP;
extern std::size_t PN;
extern std::size_t PNN;
extern std::size_t PF;
extern std::size_t PFF;
extern std::size_t PNFF;
extern std::size_t NF;
extern std::size_t NN;
extern std::size_t FF;
extern std::size_t NFF;
extern std::size_t NNFF;
extern std::size_t Fadd1;

extern Shape shape_1;
extern Shape shape_2;
extern Shape shape_X;
extern Shape shape_M;
extern Shape shape_P;
extern Shape shape_N;
extern Shape shape_F;
extern Shape shape_Fadd1;
extern Shape shape_MP;
extern Shape shape_PP;
extern Shape shape_PN;
extern Shape shape_PNN;
extern Shape shape_PF;
extern Shape shape_PFF;
extern Shape shape_PNFF;
extern Shape shape_NF;
extern Shape shape_NN;
extern Shape shape_FF;
extern Shape shape_NFF;
extern Shape shape_NNFF;

extern void static_build_shapes();
};  // namespace Dimension

namespace DATA {

namespace init {
extern VARIABLE<kids_real> Etot;
extern VARIABLE<kids_real> p;
extern VARIABLE<kids_real> x;
};  // namespace init


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
extern VARIABLE<kids_real>    sqcIA;
extern VARIABLE<kids_real>    sqcID;
extern VARIABLE<kids_real>    sqcw;
extern VARIABLE<kids_real>    sqcw0;
extern VARIABLE<kids_real>    sqcwh;

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


namespace iter {
extern VARIABLE<kids_bint> at_samplingstep_finally;
extern VARIABLE<kids_bint> at_samplingstep_initially;
extern VARIABLE<kids_real> dt;
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
};  // namespace iter


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
extern VARIABLE<kids_real> Tmod;
extern VARIABLE<kids_real> V;
extern VARIABLE<kids_int>  atoms;

namespace bath {
extern VARIABLE<kids_real> coeffs;
extern VARIABLE<kids_real> omegas;
extern VARIABLE<kids_real> p_sigma;
extern VARIABLE<kids_real> x_sigma;
};  // namespace bath


namespace coupling {
extern VARIABLE<kids_real> CL;
extern VARIABLE<kids_real> Q;
extern VARIABLE<kids_real> QL;
extern VARIABLE<kids_real> Xnj;
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

namespace rep {
extern VARIABLE<kids_real>    E;
extern VARIABLE<kids_complex> H;
extern VARIABLE<kids_real>    L;
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

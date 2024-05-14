/**@file        default_variables.h
 * @brief       provide nomination of intrinsic variables
 * @details     template of nomination of intrinsic variables, not a direct
 *              realization of variables.
 *
 * @author      Xin He
 * @date        2024-04
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
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-05-01  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef DEFAULT_VARIABLES_H
#define DEFAULT_VARIABLES_H

namespace PROJECT_NS {

// Dimension variables
namespace Dimension {
extern DSize M;  ///< No. of MonteCarlo
extern DSize P;  ///< No. of Parallel
extern DSize N;  ///< No. of Nuclear phase space pairs
extern DSize F;  ///< No. of Fock State (electronic state)

extern Shape MP({&M, &P});
extern Shape PP({&P, &P});
extern Shape PN({&P, &N});
extern Shape PNN({&P, &N, &N});
extern Shape PF({&P, &F});
extern Shape PFF({&P, &F, &F});
extern Shape PNFF({&P, &N, &F, &F});
extern Shape NF({&N, &F});
extern Shape NN({&N, &N});
extern Shape FF({&F, &F});
extern Shape NFF({&N, &F, &F});


extern int M;  ///< No. of MonteCarlo
extern int P;  ///< No. of Parallel
extern int N;  ///< No. of Nuclear phase space pairs
extern int F;  ///< No. of Fock State (electronic state)
extern int MP;
extern int PP;
extern int PN;
extern int PNN;
extern int PF;
extern int PFF;
extern int PNFF;
extern int NF;
extern int NN;
extern int FF;
extern int NFF;
extern int Fadd1;
};  // namespace Dimension

namespace DATASET {

namespace init {

extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_complex> rho_ele;
extern VARIABLE<kids_complex> rho_nuc;
extern VARIABLE<kids_real>    T;
};  // namespace init

namespace integrator {  // variables

namespace timier {
extern VARIABLE<kids_real> t;
extern VARIABLE<int>       istep;
extern VARIABLE<int>       nstep;
extern VARIABLE<int>       isamp;
extern VARIABLE<int>       nsamp;
extern VARIABLE<int>       sstep;
};  // namespace timier

extern VARIABLE<kids_real> x;
extern VARIABLE<kids_real> p;
extern VARIABLE<kids_real> m;
extern VARIABLE<kids_real> minv;
extern VARIABLE<kids_real> f;
extern VARIABLE<kids_real> fad;
extern VARIABLE<kids_real> fnad;
extern VARIABLE<kids_real> foff;
extern VARIABLE<kids_real> fadd;

extern VARIABLE<kids_complex> U;
extern VARIABLE<kids_complex> Udt;
extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_complex> c_f;
extern VARIABLE<kids_complex> c_b;
extern VARIABLE<kids_complex> rho_ele;
extern VARIABLE<kids_complex> rho_nuc;
extern VARIABLE<kids_complex> cv_mat;

extern VARIABLE<kids_real>    g;
extern VARIABLE<kids_complex> c;
extern VARIABLE<kids_real>    alpha;

extern VARIABLE<kids_complex> OpA;
extern VARIABLE<kids_complex> OpB;

extern VARIABLE<kids_real> c1;
extern VARIABLE<kids_real> c2p;

extern VARIABLE<kids_real> Ekin;
extern VARIABLE<kids_real> EPot;
extern VARIABLE<kids_real> Etot;

namespace last {
extern VARIABLE<kids_real>    x;
extern VARIABLE<kids_real>    p;
extern VARIABLE<kids_real>    grad;
extern VARIABLE<kids_real>    dV;
extern VARIABLE<kids_real>    g;
extern VARIABLE<kids_complex> c;
};  // namespace last

extern VARIABLE<int>          occ_nuc;
extern VARIABLE<kids_complex> rho_ele;
extern VARIABLE<kids_complex> rho_dual;
extern VARIABLE<kids_complex> rho_nuc;

extern VARIABLE<kids_complex> w;
extern VARIABLE<kids_complex> w_CC;
extern VARIABLE<kids_complex> w_CP;
extern VARIABLE<kids_complex> w_PP;
extern VARIABLE<kids_complex> w_AA;
extern VARIABLE<kids_complex> w_AD;
extern VARIABLE<kids_complex> w_DD;
extern VARIABLE<kids_complex> wz_A;
extern VARIABLE<kids_complex> wz_D;
extern VARIABLE<kids_complex> ww_A;
extern VARIABLE<kids_complex> ww_D;
extern VARIABLE<kids_complex> K0;
extern VARIABLE<kids_complex> K0occ;
extern VARIABLE<kids_complex> K0dia;
extern VARIABLE<kids_complex> K1;
extern VARIABLE<kids_complex> K1occ;
extern VARIABLE<kids_complex> K1dia;
extern VARIABLE<kids_complex> K2;
extern VARIABLE<kids_complex> K2occ;
extern VARIABLE<kids_complex> K2dia;
extern VARIABLE<kids_complex> K1QA;
extern VARIABLE<kids_complex> K1QAocc;
extern VARIABLE<kids_complex> K1QAdia;
extern VARIABLE<kids_complex> K2QA;
extern VARIABLE<kids_complex> K2QAocc;
extern VARIABLE<kids_complex> K2QAdia;
extern VARIABLE<kids_complex> K1DA;
extern VARIABLE<kids_complex> K1DAocc;
extern VARIABLE<kids_complex> K1DAdia;
extern VARIABLE<kids_complex> K2DA;
extern VARIABLE<kids_complex> K2DAocc;
extern VARIABLE<kids_complex> K2DAdia;
extern VARIABLE<kids_complex> K1QD;
extern VARIABLE<kids_complex> K1QDocc;
extern VARIABLE<kids_complex> K1QDdia;
extern VARIABLE<kids_complex> K2QD;
extern VARIABLE<kids_complex> K2QDocc;
extern VARIABLE<kids_complex> K2QDdia;
extern VARIABLE<kids_complex> K1DD;
extern VARIABLE<kids_complex> K1DDocc;
extern VARIABLE<kids_complex> K1DDdia;
extern VARIABLE<kids_complex> K2DD;
extern VARIABLE<kids_complex> K2DDocc;
extern VARIABLE<kids_complex> K2DDdia;
extern VARIABLE<kids_complex> OpA;
extern VARIABLE<kids_complex> OpB;
extern VARIABLE<kids_complex> TrK1A;
extern VARIABLE<kids_complex> TrK2B;

extern VARIABLE<kids_complex> U;
extern VARIABLE<kids_complex> Udt;
extern VARIABLE<kids_complex> Ubranch;

namespace tmp {
extern VARIABLE<kids_real>    direction;
extern VARIABLE<kids_real>    ve;
extern VARIABLE<kids_real>    vedE;
extern VARIABLE<kids_real>    TtTold;
extern VARIABLE<kids_real>    fproj;
extern VARIABLE<kids_real>    ftmp;
extern VARIABLE<kids_complex> wrho;
extern VARIABLE<kids_complex> invexpidiagdt;
extern VARIABLE<kids_real>    MatR_PP;
extern VARIABLE<kids_complex> MatC_PP;
extern VARIABLE<kids_complex> I_PP;
extern VARIABLE<kids_complex> fun_diag_F;
extern VARIABLE<kids_complex> fun_diag_P;
};  // namespace tmp

extern VARIABLE<kids_complex> Xcoeff;
extern VARIABLE<kids_complex> Acoeff;
extern VARIABLE<kids_complex> dtAcoeff;
extern VARIABLE<kids_complex> Hcoeff;
extern VARIABLE<kids_complex> Hbasis;
extern VARIABLE<kids_complex> UXdt;
extern VARIABLE<kids_complex> UYdt;
extern VARIABLE<kids_complex> rhored;
extern VARIABLE<kids_complex> rhored2;
extern VARIABLE<kids_complex> rhored3;
extern VARIABLE<kids_complex> Snuc;
extern VARIABLE<kids_complex> Sele;
extern VARIABLE<kids_complex> S;
extern VARIABLE<kids_complex> invS;
extern VARIABLE<kids_complex> dtlnSnuc;
extern VARIABLE<kids_complex> dtSele;

extern VARIABLE<int>       clone_count;
extern VARIABLE<int>       P_used;
extern VARIABLE<kids_real> norm;
extern VARIABLE<kids_real> veF;

namespace GWP {
extern VARIABLE<kids_real>    L;
extern VARIABLE<kids_real>    L1;
extern VARIABLE<kids_real>    L2;
extern VARIABLE<kids_complex> R;
extern VARIABLE<kids_complex> R1;
extern VARIABLE<kids_complex> R2;
extern VARIABLE<kids_complex> S1;
extern VARIABLE<kids_complex> S1h;
extern VARIABLE<kids_complex> invS1h;
extern VARIABLE<kids_complex> S2;
extern VARIABLE<kids_complex> S2h;
extern VARIABLE<kids_complex> invS2h;
extern VARIABLE<kids_complex> Sx;
};  // namespace GWP


};  // namespace integrator

namespace model {
extern VARIABLE<kids_real> mass;
extern VARIABLE<kids_real> vpes;
extern VARIABLE<kids_real> grad;
extern VARIABLE<kids_real> hess;
extern VARIABLE<kids_real> V;
extern VARIABLE<kids_real> dV;
extern VARIABLE<kids_real> ddV;
namespace rep {
extern VARIABLE<kids_real>    E;
extern VARIABLE<kids_real>    dE;
extern VARIABLE<kids_real>    T;
extern VARIABLE<kids_real>    Told;
extern VARIABLE<kids_real>    L;
extern VARIABLE<kids_complex> R;
extern VARIABLE<kids_complex> H;
};  // namespace rep
};  // namespace model

};  // namespace DATASET
};  // namespace PROJECT_NS


#endif  // DEFAULT_VARIABLES_H

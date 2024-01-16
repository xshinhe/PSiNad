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

extern VARIABLE<num_complex> c;
extern VARIABLE<num_complex> rho_ele;
extern VARIABLE<num_complex> rho_nuc;
extern VARIABLE<double> T;
};  // namespace init

namespace integrator {  // variables

namespace timier {
extern VARIABLE<double> t;
extern VARIABLE<int> istep;
extern VARIABLE<int> nstep;
extern VARIABLE<int> isamp;
extern VARIABLE<int> nsamp;
extern VARIABLE<int> sstep;
};  // namespace timier

extern VARIABLE<double> x;
extern VARIABLE<double> p;
extern VARIABLE<double> m;
extern VARIABLE<double> minv;
extern VARIABLE<double> f;
extern VARIABLE<double> fadd;
extern VARIABLE<double> g;
extern VARIABLE<num_complex> c;
extern VARIABLE<double> alpha;

extern VARIABLE<num_complex> OpA;
extern VARIABLE<num_complex> OpB;


extern VARIABLE<double> c1;
extern VARIABLE<double> c2p;

extern VARIABLE<double> Ekin;
extern VARIABLE<double> EPot;
extern VARIABLE<double> Etot;

namespace last {
extern VARIABLE<double> x;
extern VARIABLE<double> p;
extern VARIABLE<double> grad;
extern VARIABLE<double> dV;
extern VARIABLE<double> g;
extern VARIABLE<num_complex> c;
};  // namespace last

extern VARIABLE<int> occ_nuc;
extern VARIABLE<num_complex> rho_ele;
extern VARIABLE<num_complex> rho_dual;
extern VARIABLE<num_complex> rho_nuc;

extern VARIABLE<num_complex> w;
extern VARIABLE<num_complex> w_CC;
extern VARIABLE<num_complex> w_CP;
extern VARIABLE<num_complex> w_PP;
extern VARIABLE<num_complex> w_AA;
extern VARIABLE<num_complex> w_AD;
extern VARIABLE<num_complex> w_DD;
extern VARIABLE<num_complex> wz_A;
extern VARIABLE<num_complex> wz_D;
extern VARIABLE<num_complex> ww_A;
extern VARIABLE<num_complex> ww_D;
extern VARIABLE<num_complex> K0;
extern VARIABLE<num_complex> K0occ;
extern VARIABLE<num_complex> K0dia;
extern VARIABLE<num_complex> K1;
extern VARIABLE<num_complex> K1occ;
extern VARIABLE<num_complex> K1dia;
extern VARIABLE<num_complex> K2;
extern VARIABLE<num_complex> K2occ;
extern VARIABLE<num_complex> K2dia;
extern VARIABLE<num_complex> K1QA;
extern VARIABLE<num_complex> K1QAocc;
extern VARIABLE<num_complex> K1QAdia;
extern VARIABLE<num_complex> K2QA;
extern VARIABLE<num_complex> K2QAocc;
extern VARIABLE<num_complex> K2QAdia;
extern VARIABLE<num_complex> K1DA;
extern VARIABLE<num_complex> K1DAocc;
extern VARIABLE<num_complex> K1DAdia;
extern VARIABLE<num_complex> K2DA;
extern VARIABLE<num_complex> K2DAocc;
extern VARIABLE<num_complex> K2DAdia;
extern VARIABLE<num_complex> K1QD;
extern VARIABLE<num_complex> K1QDocc;
extern VARIABLE<num_complex> K1QDdia;
extern VARIABLE<num_complex> K2QD;
extern VARIABLE<num_complex> K2QDocc;
extern VARIABLE<num_complex> K2QDdia;
extern VARIABLE<num_complex> K1DD;
extern VARIABLE<num_complex> K1DDocc;
extern VARIABLE<num_complex> K1DDdia;
extern VARIABLE<num_complex> K2DD;
extern VARIABLE<num_complex> K2DDocc;
extern VARIABLE<num_complex> K2DDdia;
extern VARIABLE<num_complex> OpA;
extern VARIABLE<num_complex> OpB;
extern VARIABLE<num_complex> TrK1A;
extern VARIABLE<num_complex> TrK2B;

extern VARIABLE<num_complex> U;
extern VARIABLE<num_complex> Udt;
extern VARIABLE<num_complex> Ubranch;

namespace tmp {
extern VARIABLE<num_real> direction;
extern VARIABLE<num_real> ve;
extern VARIABLE<num_real> vedE;
extern VARIABLE<num_real> TtTold;
extern VARIABLE<num_real> fproj;
extern VARIABLE<num_real> ftmp;
extern VARIABLE<num_complex> wrho;
extern VARIABLE<num_complex> invexpidiagdt;
extern VARIABLE<num_real> MatR_PP;
extern VARIABLE<num_complex> MatC_PP;
extern VARIABLE<num_complex> I_PP;
extern VARIABLE<num_complex> fun_diag_F;
extern VARIABLE<num_complex> fun_diag_P;
};  // namespace tmp

extern VARIABLE<num_complex> Xcoeff;
extern VARIABLE<num_complex> Acoeff;
extern VARIABLE<num_complex> dtAcoeff;
extern VARIABLE<num_complex> Hcoeff;
extern VARIABLE<num_complex> Hbasis;
extern VARIABLE<num_complex> UXdt;
extern VARIABLE<num_complex> UYdt;
extern VARIABLE<num_complex> rhored;
extern VARIABLE<num_complex> rhored2;
extern VARIABLE<num_complex> rhored3;
extern VARIABLE<num_complex> Snuc;
extern VARIABLE<num_complex> Sele;
extern VARIABLE<num_complex> S;
extern VARIABLE<num_complex> invS;
extern VARIABLE<num_complex> dtlnSnuc;
extern VARIABLE<num_complex> dtSele;

extern VARIABLE<int> clone_count;
extern VARIABLE<int> P_used;
extern VARIABLE<double> norm;
extern VARIABLE<double> veF;

namespace GWP {
extern VARIABLE<double> L;
extern VARIABLE<double> L1;
extern VARIABLE<double> L2;
extern VARIABLE<num_complex> R;
extern VARIABLE<num_complex> R1;
extern VARIABLE<num_complex> R2;
extern VARIABLE<num_complex> S1;
extern VARIABLE<num_complex> S1h;
extern VARIABLE<num_complex> invS1h;
extern VARIABLE<num_complex> S2;
extern VARIABLE<num_complex> S2h;
extern VARIABLE<num_complex> invS2h;
extern VARIABLE<num_complex> Sx;
};  // namespace GWP


};  // namespace integrator

namespace model {
extern VARIABLE<double> mass;
extern VARIABLE<double> vpes;
extern VARIABLE<double> grad;
extern VARIABLE<double> hess;
extern VARIABLE<double> V;
extern VARIABLE<double> dV;
extern VARIABLE<double> ddV;
namespace rep {
extern VARIABLE<double> E;
extern VARIABLE<double> dE;
extern VARIABLE<double> T;
extern VARIABLE<double> Told;
extern VARIABLE<double> L;
extern VARIABLE<num_complex> R;
extern VARIABLE<num_complex> H;
};  // namespace rep
};  // namespace model

};  // namespace DATASET
};  // namespace PROJECT_NS


#endif  // DEFAULT_VARIABLES_H

#ifndef DEFAULT_VARIABLES_H
#define DEFAULT_VARIABLES_H

namespace kids {

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
extern VARIABLE<kids_complex> c;
extern VARIABLE<double> alpha;

extern VARIABLE<kids_complex> OpA;
extern VARIABLE<kids_complex> OpB;


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
extern VARIABLE<kids_complex> c;
};  // namespace last

extern VARIABLE<int> occ_nuc;
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
extern VARIABLE<kids_real> direction;
extern VARIABLE<kids_real> ve;
extern VARIABLE<kids_real> vedE;
extern VARIABLE<kids_real> TtTold;
extern VARIABLE<kids_real> fproj;
extern VARIABLE<kids_real> ftmp;
extern VARIABLE<kids_complex> wrho;
extern VARIABLE<kids_complex> invexpidiagdt;
extern VARIABLE<kids_real> MatR_PP;
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

extern VARIABLE<int> clone_count;
extern VARIABLE<int> P_used;
extern VARIABLE<double> norm;
extern VARIABLE<double> veF;

namespace GWP {
extern VARIABLE<double> L;
extern VARIABLE<double> L1;
extern VARIABLE<double> L2;
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
extern VARIABLE<kids_complex> R;
extern VARIABLE<kids_complex> H;
};  // namespace rep
};  // namespace model

};  // namespace DATASET
};  // namespace kids


#endif  // DEFAULT_VARIABLES_H

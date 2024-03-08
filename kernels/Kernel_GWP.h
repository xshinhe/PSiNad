#include "../core/Kernel.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

class Kernel_GWP final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_GWP"; }

    Kernel_GWP(std::shared_ptr<Kernel> kmodel, std::shared_ptr<Kernel> krepr, std::shared_ptr<Kernel> kforce)
        : _kmodel{kmodel}, _krepr{krepr}, _kforce{kforce} {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    static int calc_Ekin(kids_real* Ekin,  // [P]
                         kids_real* p,     // [P,N]
                         kids_real* m,     // [P,N]
                         int P, int N);

    static int calc_Snuc(kids_complex* S,   // [P,P]
                         kids_real* x1,     // [P,N]
                         kids_real* p1,     // [P,N]
                         kids_real* m1,     // [P,N]
                         kids_real* g1,     // [P]
                         kids_real* x2,     // [P,N]
                         kids_real* p2,     // [P,N]
                         kids_real* m2,     // [P,N]
                         kids_real* g2,     // [P]
                         kids_real* alpha,  // [N]
                         int P, int N);

    static int calc_Sele(kids_complex* S,   // [P,P]
                         kids_complex* c1,  // [P,F]
                         kids_complex* c2,  // [P,F]
                         kids_real xi,      //
                         kids_real gamma,   //
                         int P, int F);

    static int calc_dtlnSnuc(kids_complex* dtlnSnuc,  // [P,P]
                             kids_real* x,            // [P,N]
                             kids_real* p,            // [P,N]
                             kids_real* m,            // [P,N]
                             kids_real* f,            // [P,N]
                             kids_real* alpha,        // [N]
                             kids_real* Ekin,         // [P]
                             int P, int N);

    static int calc_dtSele(kids_complex* dtlnSele,  // [P,P]
                           kids_complex* Sele,      // [P,P]
                           kids_complex* c,         // [P,F]
                           kids_complex* H,         // [P,F,F]
                           kids_real* vpes,         // [P]
                           int P, int F);

    static int calc_invS(kids_complex* invS, kids_complex* S, int P);

    static double calc_density(kids_complex* rhored,  // [F,F]
                               kids_complex* Acoeff,  // [P]
                               kids_complex* Snuc,    // [P,P]
                               kids_complex* c,       // [P,F]
                               kids_real xi, kids_real gamma, int Pu, int P, int F);

    static int calc_Hbasis(kids_complex* Hbasis,  // [P,P]
                           kids_real* vpes,       // [P]
                           kids_real* grad,       // [P,N]
                           kids_real* V,          // [P,F,F]
                           kids_real* dV,         // [P,N,F,F]
                           kids_real* x,          // [P,N]
                           kids_real* p,          // [P,N]
                           kids_real* m,          // [P,N]
                           kids_real* alpha,      // [N]
                           kids_complex* Sele,    // [P,P]
                           kids_complex* c,       // [P,F]
                           int P, int N, int F    // Dimensions
    );

    static int calc_Hbasis_adia(kids_complex* Hbasis,  // [P,P]
                                kids_real* E,          // [P,F]
                                kids_real* dE,         // [P,N,F,F]
                                kids_real* x,          // [P,N]
                                kids_real* p,          // [P,N]
                                kids_real* m,          // [P,N]
                                kids_real* alpha,      // [N]
                                kids_complex* c,       // [P,F]
                                int P, int N, int F    // Dimensions
    );

   private:
    int impl_type;
    int samp_type;
    int aset_type;

    std::shared_ptr<Kernel> _kmodel;  // prepare for initial sampling
    std::shared_ptr<Kernel> _krepr;   // prepare for representation calculation
    std::shared_ptr<Kernel> _kforce;  // prepare for force calculation

    kids_real break_thres;
    int time_displace_step;
    kids_real dt;
    kids_real xi, gamma;              // for mapping kernel
    kids_real alpha0, width_scaling;  // for initial width
    kids_real *x, *p, *m, *f, *g;
    kids_real* alpha;
    kids_real* Ekin;
    kids_real* veF;

    kids_real *vpes, *grad;
    kids_real *V, *dV, *E, *dE, *T;
    kids_complex *c, *Udt, *H;

    kids_complex *Snuc, *Sele, *S, *invS;
    kids_real *L1, *L2;
    kids_complex *S1, *S1h, *invS1h, *R1;
    kids_complex *S2, *S2h, *invS2h, *R2;
    kids_complex* Sx;
    kids_complex *dtlnSnuc, *dtSele;

    ///

    kids_complex* Hbasis;
    kids_complex* Hcoeff;
    kids_complex *Acoeff, *dtAcoeff;
    kids_real* L;
    kids_complex *R, *UXdt, *UYdt, *Xcoeff;  // help for Acoeff
    kids_complex* rhored;
    kids_complex* rhored2;
    kids_complex* rhored3;

    /// temporary

    kids_real* MatR_PP;
    kids_complex* MatC_PP;
    kids_complex* I_PP;
    kids_complex* fun_diag_P;
    kids_complex* fun_diag_F;
    kids_complex* Ubranch;

    ///
    kids_real *x_last, *p_last, *grad_last, *dV_last, *g_last;
    kids_complex* c_last;

    // bool* pf_cross;
    int P_used, P_used0;
    kids_real* P_used_ptr;
    int max_clone;
    int* clone_account;
    kids_real* norm_ptr;

    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat = -1);

    int exec_kernel_impl(int stat = -1);

    int impl_0(int stat);
    int impl_1(int stat);
    int cloning();
    int death() { return 0; }
};

};  // namespace PROJECT_NS
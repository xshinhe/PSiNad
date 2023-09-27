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

    static int calc_ekin(num_real* ekin,  // [P]
                         num_real* p,     // [P,N]
                         num_real* m,     // [P,N]
                         int P, int N);

    static int calc_Snuc(num_complex* S,   // [P,P]
                         num_real* x1,     // [P,N]
                         num_real* p1,     // [P,N]
                         num_real* m1,     // [P,N]
                         num_real* g1,     // [P]
                         num_real* x2,     // [P,N]
                         num_real* p2,     // [P,N]
                         num_real* m2,     // [P,N]
                         num_real* g2,     // [P]
                         num_real* alpha,  // [N]
                         int P, int N);

    static int calc_Sele(num_complex* S,   // [P,P]
                         num_complex* c1,  // [P,F]
                         num_complex* c2,  // [P,F]
                         num_real xi,      //
                         num_real gamma,   //
                         int P, int F);

    static int calc_dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                             num_real* x,            // [P,N]
                             num_real* p,            // [P,N]
                             num_real* m,            // [P,N]
                             num_real* f,            // [P,N]
                             num_real* alpha,        // [N]
                             num_real* ekin,         // [P]
                             int P, int N);

    static int calc_dtSele(num_complex* dtlnSele,  // [P,P]
                           num_complex* Sele,      // [P,P]
                           num_complex* c,         // [P,F]
                           num_complex* H,         // [P,F,F]
                           num_real* vpes,         // [P]
                           int P, int F);

    static int calc_invS(num_complex* invS, num_complex* S, int P);

    static double calc_density(num_complex* rhored,  // [F,F]
                               num_complex* Acoeff,  // [P]
                               num_complex* Snuc,    // [P,P]
                               num_complex* c,       // [P,F]
                               num_real xi, num_real gamma, int Pu, int P, int F);

    static int calc_Hbasis(num_complex* Hbasis,  // [P,P]
                           num_real* vpes,       // [P]
                           num_real* grad,       // [P,N]
                           num_real* V,          // [P,F,F]
                           num_real* dV,         // [P,N,F,F]
                           num_real* x,          // [P,N]
                           num_real* p,          // [P,N]
                           num_real* m,          // [P,N]
                           num_real* alpha,      // [N]
                           num_complex* Sele,    // [P,P]
                           num_complex* c,       // [P,F]
                           int P, int N, int F   // Dimensions
    );

    static int calc_Hbasis_adia(num_complex* Hbasis,  // [P,P]
                                num_real* E,          // [P,F]
                                num_real* dE,         // [P,N,F,F]
                                num_real* x,          // [P,N]
                                num_real* p,          // [P,N]
                                num_real* m,          // [P,N]
                                num_real* alpha,      // [N]
                                num_complex* c,       // [P,F]
                                int P, int N, int F   // Dimensions
    );

   private:
    int impl_type;
    int samp_type;
    int aset_type;

    std::shared_ptr<Kernel> _kmodel;  // prepare for initial sampling
    std::shared_ptr<Kernel> _krepr;   // prepare for representation calculation
    std::shared_ptr<Kernel> _kforce;  // prepare for force calculation

    num_real break_thres;
    int time_displace_step;
    num_real dt;
    num_real xi, gamma;  // for mapping kernel
    num_real* gammat_ptr;
    num_real alpha0, width_scaling;  // for initial width
    num_real *x, *p, *m, *f, *g;
    num_real* alpha;
    num_real* ekin;
    num_real* veF;

    num_real *vpes, *grad;
    num_real *V, *dV, *E, *dE, *T;
    num_complex *c, *Udt, *H;

    num_complex *Snuc, *Sele, *S, *invS;
    num_real *L1, *L2;
    num_complex *S1, *S1h, *invS1h, *R1;
    num_complex *S2, *S2h, *invS2h, *R2;
    num_complex* Sx;
    num_complex *dtlnSnuc, *dtSele;

    ///

    num_complex* Hbasis;
    num_complex* Hcoeff;
    num_complex *Acoeff, *dtAcoeff;
    num_real* L;
    num_complex *R, *UXdt, *UYdt, *Xcoeff;  // help for Acoeff
    num_complex* rhored;
    num_complex* rhored2;
    num_complex* rhored3;

    /// temporary

    num_real* MatR_PP;
    num_complex* MatC_PP;
    num_complex* I_PP;
    num_complex* fun_diag_P;
    num_complex* fun_diag_F;
    num_complex* Ubranch;

    ///
    num_real *x_last, *p_last, *grad_last, *dV_last, *g_last;
    num_complex* c_last;

    // bool* pf_cross;
    int P_used, P_used0;
    num_real* P_used_ptr;
    int max_clone;
    int* clone_account;
    num_real* norm_ptr;

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
#include "../core/Kernel.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

class Kernel_GWP final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_GWP"; }

    Kernel_GWP() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    static int calc_Snuc(num_complex* S,   // [P,P]
                         num_real* x,      // [P,N]
                         num_real* p,      // [P,N]
                         num_real* alpha,  // [N]
                         int P, int N);

    static int calc_Sele(num_complex* S,  // [P,P]
                         num_complex* c,  // [P,F]
                         int P, int F);

    static int calc_dtlnSnuc(num_complex* dtlnSnuc,  // [P,P]
                             num_real* x,            // [P,N]
                             num_real* p,            // [P,N]
                             num_real* m,            // [P,N]
                             num_real* f,            // [P,N]
                             num_real* alpha,        // [N]
                             int P, int N);

    static int calc_dtlnSele(num_complex* dtlnSele,  // [P,P]
                             num_complex* c,         // [P,F]
                             num_complex* dUdt,      // [P,F,F]
                             int P, int F);

    static int calc_invS(num_complex* invS, num_complex* S, int P);

    static double calc_density(num_complex* rhored,  // [F,F]
                               num_complex* Acoeff,  // [P]
                               num_complex* Snuc,    // [P,P]
                               num_complex* c,       // [P,F]
                               int P, int F);

    static int calc_Hbasis(num_complex* Hbasis,  // [P,P]
                           num_real* V,          // [P,F,F]
                           num_real* dV,         // [P,N,F,F]
                           num_real* x,          // [P,N]
                           num_real* p,          // [P,N]
                           num_real* m,          // [P,N]
                           num_real* alpha,      // [N]
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
    num_real dt;
    num_real xi, gamma, alpha0;
    num_real *alpha, *x, *p, *m, *f;
    num_real *V, *dV, *E, *dE;
    num_complex *c, *dUdt;

    num_complex *S, *invS;
    num_complex *Snuc, *Sele;
    num_complex *dtlnSnuc, *dtlnSele;
    num_complex *S1, *S2;

    num_complex* Hbasis;
    num_complex* Hcoeff;
    num_complex *Acoeff, *dtAcoeff;
    num_complex* rhored;

    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat = -1);

    int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS
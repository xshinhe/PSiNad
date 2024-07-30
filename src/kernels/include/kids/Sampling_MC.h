#include "kids/Kernel.h"

namespace PROJECT_NS {

class Sampling_MC final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Sampling_MC(std::shared_ptr<Kernel> kmodel, std::shared_ptr<Kernel> krepr, std::shared_ptr<Kernel> kforce)
        : _kmodel{kmodel}, _krepr{krepr}, _kforce{kforce} {};

   private:
    int impl_type;
    int samp_type;
    int aset_type;

    std::shared_ptr<Kernel> _kmodel;  // prepare for initial sampling
    std::shared_ptr<Kernel> _krepr;   // prepare for representation calculation
    std::shared_ptr<Kernel> _kforce;  // prepare for force calculation

    kids_real  break_thres;
    int        time_displace_step;
    kids_real  dt;
    kids_real  xi, gamma;              // for mapping kernel
    kids_real  alpha0, width_scaling;  // for initial width
    kids_real *x, *p, *m, *f, *g;
    kids_real* alpha;
    kids_real* Ekin;
    kids_real *ve, *veF;

    kids_real *   vpes, *grad;
    kids_real *   V, *dV, *eig, *dE, *T;
    kids_complex *c, *U, *Udt, *H;

    kids_complex *Snuc, *Sele, *S, *invS;
    kids_real *   L1, *L2;
    kids_complex *S1, *S1h, *invS1h, *R1;
    kids_complex *S2, *S2h, *invS2h, *R2;
    kids_complex* Sx;
    kids_complex *dtlnSnuc, *dtSele;

    ///

    kids_complex* Hbasis;
    kids_complex* Hcoeff;
    kids_complex *Acoeff, *dtAcoeff;
    kids_real*    L;
    kids_complex *R, *UXdt, *UYdt, *Xcoeff;  // help for Acoeff
    kids_complex* rhored;
    kids_complex* rhored2;
    kids_complex* rhored3;

    /// temporary

    kids_real*    MatR_PP;
    kids_complex* MatC_PP;
    kids_complex* I_PP;
    kids_complex* fun_diag_P;
    kids_complex* fun_diag_F;
    kids_complex* Ubranch;

    ///
    kids_real *   x_last, *p_last, *grad_last, *dV_last, *g_last;
    kids_complex *c_last, *c_init;

    // bool* pf_cross;
    int        P_used, P_used0;
    kids_real* P_used_ptr;
    int        max_clone;
    int*       clone_account;
    kids_real* norm_ptr;

    kids_complex* w;
    kids_complex *rho_nuc, *rho_ele;
    kids_int*     occ_nuc;
    kids_real*    T_init;

    int occ0;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    // Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS
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

    kids_real       break_thres;
    int             time_displace_step;
    kids_real       dt;
    kids_real       xi, gamma;              // for mapping kernel
    kids_real       alpha0, width_scaling;  // for initial width
    span<kids_real> x, p, m, f, g;
    span<kids_real> alpha;
    span<kids_real> Ekin;
    span<kids_real> ve, veF;

    span<kids_real>    vpes, grad;
    span<kids_real>    V, dV, eig, dE, T;
    span<kids_complex> c, U, Udt, H;

    span<kids_complex> Snuc, Sele, S, invS;
    span<kids_real>    L1, L2;
    span<kids_complex> S1, S1h, invS1h, R1;
    span<kids_complex> S2, S2h, invS2h, R2;
    span<kids_complex> Sx;
    span<kids_complex> dtlnSnuc, dtSele;

    ///

    span<kids_complex> Hbasis;
    span<kids_complex> Hcoeff;
    span<kids_complex> Acoeff, dtAcoeff;
    span<kids_real>    L;
    span<kids_complex> R, UXdt, UYdt, Xcoeff;  // help for Acoeff
    span<kids_complex> rhored;
    span<kids_complex> rhored2;
    span<kids_complex> rhored3;

    /// temporary

    span<kids_real>    MatR_PP;
    span<kids_complex> MatC_PP;
    span<kids_complex> I_PP;
    span<kids_complex> fun_diag_P;
    span<kids_complex> fun_diag_F;
    span<kids_complex> Ubranch;

    ///
    span<kids_real>    x_last, p_last, grad_last, dV_last, g_last;
    span<kids_complex> c_last, c_init;

    // bool* pf_cross;
    int             P_used, P_used0;
    span<kids_int>  P_used_ptr;
    int             max_clone;
    span<kids_int>  clone_account;
    span<kids_real> norm_ptr;

    span<kids_complex> w;
    span<kids_complex> rho_nuc, rho_ele;
    span<kids_int>     occ_nuc;
    span<kids_real>    T_init;

    int occ0;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    // Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS
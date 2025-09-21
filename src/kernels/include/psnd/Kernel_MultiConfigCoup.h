#include "psnd/Kernel.h"

namespace PROJECT_NS {

class Kernel_MultiConfigCoup final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_MultiConfigCoup(std::shared_ptr<Kernel> kmodel, std::shared_ptr<Kernel> krepr,
                           std::shared_ptr<Kernel> kforce)
        : _kmodel{kmodel}, _krepr{krepr}, _kforce{kforce} {};

    static int calc_Ekin(span<psnd_real> Ekin,  // [P]
                         span<psnd_real> p,     // [P,N]
                         span<psnd_real> m,     // [P,N]
                         int P, int N);

    static int calc_Snuc(span<psnd_complex> S,      // [P,P]
                         span<psnd_real>    x1,     // [P,N]
                         span<psnd_real>    p1,     // [P,N]
                         span<psnd_real>    m1,     // [P,N]
                         span<psnd_real>    g1,     // [P]
                         span<psnd_real>    x2,     // [P,N]
                         span<psnd_real>    p2,     // [P,N]
                         span<psnd_real>    m2,     // [P,N]
                         span<psnd_real>    g2,     // [P]
                         span<psnd_real>    alpha,  // [N]
                         int P, int N);

    static int calc_Sele(span<psnd_complex> S,      // [P,P]
                         span<psnd_complex> c1,     // [P,F]
                         span<psnd_complex> c2,     // [P,F]
                         psnd_real          xi,     //
                         psnd_real          gamma,  //
                         int P, int F);

    static int calc_dtlnSnuc(span<psnd_complex> dtlnSnuc,  // [P,P]
                             span<psnd_real>    x,         // [P,N]
                             span<psnd_real>    p,         // [P,N]
                             span<psnd_real>    m,         // [P,N]
                             span<psnd_real>    f,         // [P,N]
                             span<psnd_real>    alpha,     // [N]
                             span<psnd_real>    Ekin,      // [P]
                             int P, int N);

    static int calc_dtSele(span<psnd_complex> dtlnSele,  // [P,P]
                           span<psnd_complex> Sele,      // [P,P]
                           span<psnd_complex> c,         // [P,F]
                           span<psnd_complex> H,         // [P,F,F]
                           span<psnd_real>    vpes,      // [P]
                           int P, int F);

    static int calc_invS(span<psnd_complex> invS, span<psnd_complex> S, int P);

    static double calc_density(span<psnd_complex> rhored,  // [F,F]
                               span<psnd_complex> Acoeff,  // [P]
                               span<psnd_complex> Snuc,    // [P,P]
                               span<psnd_complex> c,       // [P,F]
                               span<psnd_complex> Mtmp,    // [P,P]
                               psnd_real xi, psnd_real gamma, int Pu, int P, int F);

    static int calc_Hbasis(span<psnd_complex> Hbasis,  // [P,P]
                           span<psnd_real>    vpes,    // [P]
                           span<psnd_real>    grad,    // [P,N]
                           span<psnd_real>    V,       // [P,F,F]
                           span<psnd_real>    dV,      // [P,N,F,F]
                           span<psnd_real>    x,       // [P,N]
                           span<psnd_real>    p,       // [P,N]
                           span<psnd_real>    m,       // [P,N]
                           span<psnd_real>    alpha,   // [N]
                           span<psnd_complex> Sele,    // [P,P]
                           span<psnd_complex> c,       // [P,F]
                           int P, int N, int F         // Dimensions
    );

    static int calc_Hbasis_adia(span<psnd_complex> Hbasis,  // [P,P]
                                span<psnd_real>    E,       // [P,F]
                                span<psnd_real>    dE,      // [P,N,F,F]
                                span<psnd_real>    x,       // [P,N]
                                span<psnd_real>    p,       // [P,N]
                                span<psnd_real>    m,       // [P,N]
                                span<psnd_real>    alpha,   // [N]
                                span<psnd_complex> c,       // [P,F]
                                int P, int N, int F         // Dimensions
    );

   private:
    int impl_type;
    int samp_type;
    int aset_type;

    std::shared_ptr<Kernel> _kmodel;  // prepare for initial sampling
    std::shared_ptr<Kernel> _krepr;   // prepare for representation calculation
    std::shared_ptr<Kernel> _kforce;  // prepare for force calculation

    psnd_real       break_thres;
    int             time_displace_step;
    psnd_real       dt;
    psnd_real       xi, gamma;              // for mapping kernel
    psnd_real       alpha0, width_scaling;  // for initial width
    span<psnd_real> x, p, m, f, g;
    span<psnd_real> alpha;
    span<psnd_real> Ekin;
    span<psnd_real> ve, veF;

    span<psnd_real>    vpes, grad;
    span<psnd_real>    V, dV, eig, dE, T;
    span<psnd_complex> c, U, Udt, H;

    span<psnd_complex> Snuc, Sele, S, invS;
    span<psnd_real>    L1, L2;
    span<psnd_complex> S1, S1h, invS1h, R1;
    span<psnd_complex> S2, S2h, invS2h, R2;
    span<psnd_complex> Sx;
    span<psnd_complex> dtlnSnuc, dtSele;

    ///

    span<psnd_complex> Hbasis;
    span<psnd_complex> Hcoeff;
    span<psnd_complex> Acoeff, dtAcoeff;
    span<psnd_real>    L;
    span<psnd_complex> R, UXdt, UYdt, Xcoeff;  // help for Acoeff
    span<psnd_complex> rhored;
    span<psnd_complex> rhored2;
    span<psnd_complex> rhored3;

    /// temporary

    span<psnd_real>    MatR_PP;
    span<psnd_complex> MatC_PP;
    span<psnd_complex> I_PP;
    span<psnd_complex> fun_diag_P;
    span<psnd_complex> fun_diag_F;
    span<psnd_complex> Ubranch;

    ///
    span<psnd_real>    x_last, p_last, grad_last, dV_last, g_last;
    span<psnd_complex> c_last, c_init;

    // bool* pf_cross;
    int             P_used, P_used0;
    span<psnd_int>  P_used_ptr;
    int             max_clone;
    span<psnd_int>  clone_account;
    span<psnd_real> norm_ptr;

    span<psnd_complex> w;
    span<psnd_complex> rho_nuc;
    span<psnd_int>     occ_nuc;
    span<psnd_real>    T_init;

    int occ0;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);

    Status& impl_0(Status& stat);
    Status& impl_1(Status& stat);
    int     cloning();
    int     death() { return 0; }
};

};  // namespace PROJECT_NS
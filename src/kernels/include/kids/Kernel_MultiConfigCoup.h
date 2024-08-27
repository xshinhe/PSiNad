#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_MultiConfigCoup final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_MultiConfigCoup(std::shared_ptr<Kernel> kmodel, std::shared_ptr<Kernel> krepr,
                           std::shared_ptr<Kernel> kforce)
        : _kmodel{kmodel}, _krepr{krepr}, _kforce{kforce} {};

    static int calc_Ekin(span<kids_real> Ekin,  // [P]
                         span<kids_real> p,     // [P,N]
                         span<kids_real> m,     // [P,N]
                         int P, int N);

    static int calc_Snuc(span<kids_complex> S,      // [P,P]
                         span<kids_real>    x1,     // [P,N]
                         span<kids_real>    p1,     // [P,N]
                         span<kids_real>    m1,     // [P,N]
                         span<kids_real>    g1,     // [P]
                         span<kids_real>    x2,     // [P,N]
                         span<kids_real>    p2,     // [P,N]
                         span<kids_real>    m2,     // [P,N]
                         span<kids_real>    g2,     // [P]
                         span<kids_real>    alpha,  // [N]
                         int P, int N);

    static int calc_Sele(span<kids_complex> S,      // [P,P]
                         span<kids_complex> c1,     // [P,F]
                         span<kids_complex> c2,     // [P,F]
                         kids_real          xi,     //
                         kids_real          gamma,  //
                         int P, int F);

    static int calc_dtlnSnuc(span<kids_complex> dtlnSnuc,  // [P,P]
                             span<kids_real>    x,         // [P,N]
                             span<kids_real>    p,         // [P,N]
                             span<kids_real>    m,         // [P,N]
                             span<kids_real>    f,         // [P,N]
                             span<kids_real>    alpha,     // [N]
                             span<kids_real>    Ekin,      // [P]
                             int P, int N);

    static int calc_dtSele(span<kids_complex> dtlnSele,  // [P,P]
                           span<kids_complex> Sele,      // [P,P]
                           span<kids_complex> c,         // [P,F]
                           span<kids_complex> H,         // [P,F,F]
                           span<kids_real>    vpes,      // [P]
                           int P, int F);

    static int calc_invS(span<kids_complex> invS, span<kids_complex> S, int P);

    static double calc_density(span<kids_complex> rhored,  // [F,F]
                               span<kids_complex> Acoeff,  // [P]
                               span<kids_complex> Snuc,    // [P,P]
                               span<kids_complex> c,       // [P,F]
                               span<kids_complex> Mtmp,    // [P,P]
                               kids_real xi, kids_real gamma, int Pu, int P, int F);

    static int calc_Hbasis(span<kids_complex> Hbasis,  // [P,P]
                           span<kids_real>    vpes,    // [P]
                           span<kids_real>    grad,    // [P,N]
                           span<kids_real>    V,       // [P,F,F]
                           span<kids_real>    dV,      // [P,N,F,F]
                           span<kids_real>    x,       // [P,N]
                           span<kids_real>    p,       // [P,N]
                           span<kids_real>    m,       // [P,N]
                           span<kids_real>    alpha,   // [N]
                           span<kids_complex> Sele,    // [P,P]
                           span<kids_complex> c,       // [P,F]
                           int P, int N, int F         // Dimensions
    );

    static int calc_Hbasis_adia(span<kids_complex> Hbasis,  // [P,P]
                                span<kids_real>    E,       // [P,F]
                                span<kids_real>    dE,      // [P,N,F,F]
                                span<kids_real>    x,       // [P,N]
                                span<kids_real>    p,       // [P,N]
                                span<kids_real>    m,       // [P,N]
                                span<kids_real>    alpha,   // [N]
                                span<kids_complex> c,       // [P,F]
                                int P, int N, int F         // Dimensions
    );

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
    span<kids_complex> rho_nuc;
    span<kids_int>     occ_nuc;
    span<kids_real>    T_init;

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
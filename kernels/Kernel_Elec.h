#ifndef Kernel_Elec_H
#define Kernel_Elec_H

#include "../core/Kernel.h"

namespace PROJECT_NS {
/**
 * @brief Kernel_Elec:
 *    1) evolve by rho_ele
 *    2) evolve by c(:,r), rank-decomposition of density rho_ele
 */
class Kernel_Elec final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec"; }

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_c(kids_complex *ker, kids_complex *c, kids_real xi, kids_real gamma, int fdim);

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_rho(kids_complex *ker, kids_complex *rho, kids_real xi, kids_real gamma, int fdim,
                            bool quantize = false, int occ = -1);

    /**
     * read parameters
     */
    static int occ0;

    /**
     * dynamics variables for electronic DOFs
     * @note: the evolution of U refers "Kernel_Update.h"
     *   c = U*c_init
     *   rho_ele = U*rho_ele*U^
     */
    static kids_complex *U;                         ///< propogator along classical path
    static kids_complex *c, *c_init;                ///< electronic vector
    static kids_complex *rho_ele, *rho_ele_init;    ///< electronic density
    static kids_complex *rho_dual, *rho_dual_init;  ///< electronic density
    static kids_real *T, *T_init;

    /**
     * weighting density for nuclear force
     */
    static int *occ_nuc;
    static kids_complex *rho_nuc, *rho_nuc_init;

    /**
     * kernels for time correlation function
     */
    static kids_complex *w;  ///< initial measurement of the phase point
    static kids_complex *wz_A;
    static kids_complex *wz_D;
    static kids_complex *ww_A;
    static kids_complex *ww_D;
    static kids_complex *ww_A_init;
    static kids_complex *ww_D_init;
    static kids_complex *w_AA, *w_AD, *w_DD, *w_CC, *w_CP, *w_PP;
    static kids_complex *K0, *K0occ, *K0dia;        ///< partial version of K0
    static kids_complex *K1, *K1occ, *K1dia;        ///< partial version of K1
    static kids_complex *K2, *K2occ, *K2dia;        ///< partial version of K2
    static kids_complex *K1QA, *K1QAocc, *K1QAdia;  ///< Simplex Quantization
    static kids_complex *K2QA, *K2QAocc, *K2QAdia;  ///< Heaviside Quantization
    static kids_complex *K1DA, *K1DAocc, *K1DAdia;
    static kids_complex *K2DA, *K2DAocc, *K2DAdia;
    static kids_complex *K1QD, *K1QDocc, *K1QDdia;  ///< Simplex Quantization
    static kids_complex *K2QD, *K2QDocc, *K2QDdia;  ///< Heaviside Quantization
    static kids_complex *K1DD, *K1DDocc, *K1DDdia;
    static kids_complex *K2DD, *K2DDocc, *K2DDdia;

    static kids_complex *OpA, *OpB;
    static kids_complex *TrK1A, *TrK2B;

   private:
    void read_param_impl(Param *PM);

    void init_data_impl(DataSet *DS);

    void init_calc_impl(int stat = -1);

    int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_H

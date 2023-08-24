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
    static int ker_from_c(num_complex *ker, num_complex *c, num_real xi, num_real gamma, int fdim);

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_rho(num_complex *ker, num_complex *rho, num_real xi, num_real gamma, int fdim,
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
    static num_complex *U;                       ///< propogator along classical path
    static num_complex *c, *c_init;              ///< electronic vector
    static num_complex *rho_ele, *rho_ele_init;  ///< electronic density

    /**
     * weighting density for nuclear force
     */
    static int *occ_nuc;
    static num_complex *rho_nuc, *rho_nuc_init;

    /**
     * kernels for time correlation function
     */
    static num_complex *w;                   ///< initial measurement of the phase point
    static num_complex *K1, *K1occ, *K1dia;  ///< partial version of K1
    static num_complex *K2, *K2occ, *K2dia;  ///< partial version of K2
    static num_complex *K1D, *K1Docc, *K1Ddia;
    static num_complex *K2D, *K2Docc, *K2Ddia;

    static num_complex *OpA, *OpB;
    static num_complex *TrK1A, *TrK2B;

   private:
    void read_param_impl(Param *PM);

    void init_data_impl(DataSet *DS);

    void init_calc_impl(int stat = -1);

    int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_H

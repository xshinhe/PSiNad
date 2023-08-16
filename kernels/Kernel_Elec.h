#ifndef Kernel_Elec_H
#define Kernel_Elec_H

#include "../core/Kernel.h"

namespace PROJECT_NS {
/**
 * @brief Kernel_Elec : An interface for transforming electonic DOFs
 */
class Kernel_Elec final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec"; }

    /**
     * @brief convert mapping variables to c (electonic amplititude)
     */
    static int c_from_mapvar(num_complex* c, num_real* mapvar, int fdim);
    /**
     * @brief convert c (electonic amplititude) to mapping variables
     */
    static int mapvar_from_c(num_real* mapvar, num_complex* c, int fdim);

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_c(num_complex* ker, num_complex* c, num_real xi, num_real gamma, int fdim);

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_rho(num_complex* ker, num_complex* rho, num_real xi, num_real gamma, int fdim);

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_rho_quantize(num_complex* ker, num_complex* rho, num_real xi, num_real gamma, int occt,
                                     int fdim);

    /**
     * @brief sampling mapping variables from focused condition (A & Q as inputs)
     */
    static int mapvar_focus(num_real* mapvar, int fdim);

    static num_complex* w;  // measure of phase point
    static num_complex* c;
    static num_real* mapvar;
    static num_real* mapx;
    static num_real* mapp;
    static num_real* mapA;
    static num_real* mapQ;
    static int occ0;

    // two densities for dynamic
    static int* occ_nuc;
    static num_complex* rho_ele;  // electronic density
    static num_complex* gmat;
    static num_complex* rho_nuc;  // nuclear weighting density

    // time correlation function
    static num_complex *K1, *wK1, *wK1occ, *wK1dia;
    static num_complex *K2, *K2dia;
    static num_complex* K1Q;  // quantized K1, or alternative K1
    static num_complex* K2Q;  // quantized K2, or alternative K2

    static num_complex *OpA, *OpB, *TrK1A, *TrK2B;
    // friend class Kernel_Elec_CMM;
   private:
    void read_param_impl(Param* PM);
    void init_data_impl(DataSet* DS);
    void init_calc_impl(int stat = -1);
    int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_H

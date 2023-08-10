/**
 * @file Kernel_Elec.h
 * @author xshinhe
 * @version 1.1
 * @date 2023-03
 * @brief initialization kernels for electonic DOFs
 * @details
 *  The initialization of electonic DOFs are tightly related to Solver.
 *  Use it in Solver's Kernel_Builder();
 */

#ifndef Kernel_Elec_LSC_H
#define Kernel_Elec_LSC_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in LSC
 */
class Kernel_Elec_LSC final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_LSC"; }

    /**
     * @brief function for gamma_wigner: gamma = (sqrt(1+F)-1)/F
     */
    static double gamma_wigner(int fdim);

    /**
     * @brief sampling mapping variables from gaussian distribution
     */
    static int mapvar_gauss(num_real *mapvar, num_real variance, int fdim) {
        Kernel_Random::rand_gaussian(mapvar, 2 * fdim, sqrt(variance));
        return 0;
    }

   private:
    num_complex *w;         // measure of phase point mu(dX) = w(X)dX
    num_complex *c, *gmat;  // only used for check
    num_complex *rho, *rho_nuc;
    num_complex *K0, *wK0, *wK0occ, *wK0dia;
    num_complex *Kt, *Ktdia;

    num_real *mapvar, *mapx, *mapp, *mapA, *mapQ;

    int occ;
    num_real gamma0, gammat, xi0, xit, totact;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // Kernel_Elec_LSC_H

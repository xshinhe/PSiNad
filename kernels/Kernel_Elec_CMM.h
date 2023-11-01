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

#ifndef Kernel_Elec_CMM_H
#define Kernel_Elec_CMM_H

#include "../core/Kernel.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in CMM
 */
class Kernel_Elec_CMM final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_CMM"; }

    Kernel_Elec_CMM() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    /**
     * @brief function for: gamma_wigner(F) = (sqrt(1+F)-1)/F
     */
    static double gamma_wigner(int fdim);

    /**
     * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
     */
    static int c_sphere(num_complex *c, int fdim);

    static int c_focus(num_complex *c, double xi, double gamma, int occ, int fdim);

   private:
    num_real gamma1, gamma2, xi1, xi2;
    bool use_cv  = false;
    bool use_wmm = false;  // in this case, gamma1 will be used as delta in wMM

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // Kernel_Elec_CMM_H

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

#ifndef Kernel_Elec_SH_H
#define Kernel_Elec_SH_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in SH
 */
class Kernel_Elec_SH final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_SH"; }

    Kernel_Elec_SH() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    /**
     * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
     */
    static int mapvar_sphere(num_real *mapvar, num_real Rc2, int fdim);

   private:
    int occ;
    num_real gamma0, gammat, xi0, xit, totact;

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // Kernel_Elec_SH_H

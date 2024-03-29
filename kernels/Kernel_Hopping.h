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

#ifndef Kernel_Hopping_H
#define Kernel_Hopping_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in SH
 */
class Kernel_Hopping final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Hopping"; }

    Kernel_Hopping() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    // @brief: generate hopping state from iocc (but with change current state)
    static int max_choose(kids_complex* rho);

    static int pop_choose(kids_complex* rho);

    static int pop_neg_choose(kids_complex* rho);

    // @brief: generate hopping state from iocc (but with change current state)
    static int hopping_choose(kids_complex* rho, kids_complex* H, int from, kids_real dt);

    static void hopping_direction(kids_real* direction, kids_real* dE, int from, int to);

    static int hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm, kids_real* E,  //
                               int from, int to, bool reflect);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Hopping_H

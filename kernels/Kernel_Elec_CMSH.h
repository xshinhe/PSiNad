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

#ifndef Kernel_Elec_CMSH_H
#define Kernel_Elec_CMSH_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in CMSH
 */
class Kernel_Elec_CMSH final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_CMSH"; }

    Kernel_Elec_CMSH() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

   private:
    num_real gamma1, gamma2, xi1, xi2;
    num_real alpha;
    bool use_cv         = true;
    bool use_wmm        = false;  // in this case, gamma1 will be used as delta in wMM
    bool reflect        = true;
    bool conserve_scale = false;
    int hopping_type1;
    int hopping_type2;
    int hopping_type3;

    double dt;
    num_real* p;
    num_real* m;
    num_real* vpes;
    num_real* Epot;
    num_real* Ekin;
    num_real* Etot;
    num_real* Etot_init;
    num_real* fadd;
    num_real* direction;
    num_real *E, *dE, *T;
    num_complex* H;


    int hopping_impulse(num_real* direction, num_real* np, num_real* nm,  //
                        num_real Efrom, num_real Eto, int from, int to, bool reflect);

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_CMSH_H

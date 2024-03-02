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

#ifndef Kernel_Elec_MMSH_H
#define Kernel_Elec_MMSH_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace kids {

DEFINE_POLICY(MMSHPolicy,
              MASH1,  // arXiv:2212.11773
              MASH2,  // arXiv:2305.08835
              MFSH1);


/**
 * @brief initialization kernel for electonic DOFs in MMSH
 */
class Kernel_Elec_MMSH final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_MMSH"; }

    Kernel_Elec_MMSH() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    static inline double gamma_opt(int fdim) {
        double sum = 0.0f;
        for (int i = 0; i < fdim; ++i) sum += 1.0f / (i + 1);
        return ((fdim - sum) / (sum - 1.0f)) / fdim;
    }

    static void hopping_direction(kids_real* direction, kids_real* E, kids_real* dE, kids_complex* rho, int from,
                                  int to);

   private:
    MMSHPolicy::_type mmsh_type;
    bool sumover;
    bool focused;
    bool reflect;
    bool hopping;

    double xi;
    double gamma;
    bool use_cv = false;


    double dt;
    kids_real* p;
    kids_real* m;
    kids_real* direction;
    kids_real *E, *dE, *T;
    kids_complex* H;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids

#endif  // Kernel_Elec_MMSH_H

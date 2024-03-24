#ifndef Kernel_Iter_H
#define Kernel_Iter_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * @brief iterative kernel wrapper/(interface) for other kernels
 */
class Kernel_Iter final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter"; }

   private:
    double t0, *t0_ptr;
    double t, *t_ptr;
    double dt, *dt_ptr;
    double tend, *tend_ptr;
    double tsec, *tsec_ptr;
    int *succ_ptr;
    bool *at_samplingstep_initially_ptr;
    bool *at_samplingstep_finally_ptr;

    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_H

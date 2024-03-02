#ifndef Kernel_Timer_H
#define Kernel_Timer_H

#include "../core/Kernel.h"

namespace kids {

class Kernel_Timer final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Timer"; }

   private:
    double t, *t_ptr;
    double t0, tend, dt;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids


#endif  // Kernel_Timer_H

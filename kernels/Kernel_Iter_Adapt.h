#ifndef Kernel_Iter_Adapt_H
#define Kernel_Iter_Adapt_H

#include "../core/Kernel.h"

namespace kids {

class Kernel_Iter_Adapt final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter_Adapt"; }

   private:
    // KData<OnlyRead> t0; // read from Param
    // KData<OnlyRead> tend;
    // KData<Register> dt;
    // KData<OnlyRefs> other;

    double t0, *t0_ptr;
    double t, *t_ptr;
    double dt, *dt_ptr;
    double tend, *tend_ptr;
    double tsec, *tsec_ptr;
    int *succ_ptr;
    int *do_recd_ptr;
    int *do_prec_ptr;

    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

    const std::vector<std::string> backup_fields = {"x", "p", "rho_ele", "rho_nuc", "U", "occ"};

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual void init_data_impl(DataSet *DS);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids


#endif  // Kernel_Iter_Adapt_H

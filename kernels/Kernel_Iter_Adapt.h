#ifndef Kernel_Iter_Adapt_H
#define Kernel_Iter_Adapt_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

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
    bool *succ_ptr;
    bool *frez_ptr;
    bool *do_recd_ptr;
    bool *do_prec_ptr;

    int msize, *tsize_ptr, *dtsize_ptr;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;
    int nbackup;

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin"};

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual void init_data_impl(DataSet *DS);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_Adapt_H

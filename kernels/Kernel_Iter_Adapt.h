#ifndef Kernel_Iter_Adapt_H
#define Kernel_Iter_Adapt_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Iter_Adapt final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter_Adapt"; }

   private:
    double t0, *t0_ptr;
    double t, *t_ptr;
    double dt, *dt_ptr;
    double tend, *tend_ptr;
    bool* at_samplingstep_initially_ptr;
    bool* at_samplingstep_finally_ptr;

    bool* succ_ptr;
    bool* frez_ptr;
    bool* last_attempt_ptr;
    int* fail_type_ptr;  // record the failure information (longtime keeped)

    int msize, *tsize_ptr, *dtsize_ptr;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;
    int nbackup;

    double time_unit;

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin", "Epot"};

    virtual void read_param_impl(Param* PM);

    virtual void init_calc_impl(int stat = -1);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_Adapt_H

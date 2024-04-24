#ifndef Kernel_Iter_Adapt_H
#define Kernel_Iter_Adapt_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Iter_Adapt final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter_Adapt"; }

   private:
    double     t0, *t0_ptr;
    double     t, *t_ptr;
    double     dt, *dt_ptr;
    double     tend, *tend_ptr;
    kids_bint* at_samplingstep_initially_ptr;
    kids_bint* at_samplingstep_finally_ptr;

    kids_bint* succ_ptr;
    kids_bint* frez_ptr;
    kids_bint* last_attempt_ptr;
    int*       fail_type_ptr;  // record the failure information (longtime keeped)

    int msize, *tsize_ptr, *dtsize_ptr;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;
    int nbackup;

    double time_unit;

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin", "Epot"};

    virtual void setInputParam_impl(std::shared_ptr<Param>& PM);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_Adapt_H

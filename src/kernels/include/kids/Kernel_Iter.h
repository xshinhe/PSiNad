#ifndef Kernel_Iter_H
#define Kernel_Iter_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * Minimal iterator for integration of the equations of motion
 */
class Kernel_Iter final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    double     t0, *t0_ptr;
    double     t, *t_ptr;
    double     dt, *dt_ptr;
    double     tend, *tend_ptr;
    double     tsec, *tsec_ptr;
    int*       succ_ptr;
    kids_bint* at_samplingstep_initially_ptr;
    kids_bint* at_samplingstep_finally_ptr;

    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_H

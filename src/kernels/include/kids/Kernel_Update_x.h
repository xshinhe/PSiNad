#ifndef Kernel_Update_x_H
#define Kernel_Update_x_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_x final : public Kernel {
   public:
    Kernel_Update_x(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_x"; }

   private:
    double *x, *p, *m, *minv;

    double     scale, *dt_ptr;
    kids_bint* frez_ptr;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Update_x_H

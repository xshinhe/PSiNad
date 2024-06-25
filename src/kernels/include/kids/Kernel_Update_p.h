#ifndef Kernel_Update_p_H
#define Kernel_Update_p_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_p : public Kernel {
   public:
    Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

   private:
    double *   p, *f, *minv;
    double*    Ekin;
    double     scale, *dt_ptr;
    kids_bint* frez_ptr;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_p_H

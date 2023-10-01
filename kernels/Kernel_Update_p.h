#ifndef Kernel_Update_p_H
#define Kernel_Update_p_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_p : public Kernel {
   public:
    Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_p"; }

   private:
    double *p, *f;

    double scale;
    double dt, sdt;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_p_H

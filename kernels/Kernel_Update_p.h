#ifndef Kernel_Update_p_H
#define Kernel_Update_p_H

#include "../core/Kernel.h"

namespace kids {

class Kernel_Update_p : public Kernel {
   public:
    Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_p"; }

   private:
    double *p, *f, *minv;
    double* Ekin;
    double scale, *dt_ptr;

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace kids


#endif  // Kernel_Update_p_H

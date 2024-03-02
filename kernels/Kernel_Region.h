#ifndef Kernel_Region_H
#define Kernel_Region_H

#include "../core/Kernel.h"
#include "Kernel_Elec.h"

namespace kids {

class Kernel_Region final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Region"; }

   private:
    region_type;

    kids_real *f, *grad, *dV, *dE, *Force, *T;

    virtual void read_param_impl(Param* PM){};

    virtual void init_data_impl(DataSet* DS){};

    virtual void init_calc_impl(int stat = -1){};

    virtual int exec_kernel_impl(int stat = -1) { return 0; }
};

};  // namespace kids

#endif  // Kernel_Region_H
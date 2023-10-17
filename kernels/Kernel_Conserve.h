#ifndef Kernel_Conserve_H
#define Kernel_Conserve_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Conserve final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Conserve"; }

   private:
    num_real* E;
    num_real* p;
    num_real* m;
    num_real* Etot;
    num_real* Etot_init;
    num_real* Ekin;
    num_real* Epot;
    num_real* vpes;
    double conserve_scale;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Conserve_H

#ifndef Kernel_Conserve_H
#define Kernel_Conserve_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Conserve final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Conserve"; }

   private:
    kids_real* E;
    kids_real* p;
    kids_real* m;
    kids_real* Etot;
    kids_real* Etot_prev;
    kids_real* Etot_init;
    kids_real* Ekin;
    kids_real* Epot;
    kids_real* vpes;

    bool conserve_scale;
    bool* succ_ptr;
    bool* frez_ptr;
    bool* last_attempt_ptr;
    int* fail_type_ptr;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Conserve_H

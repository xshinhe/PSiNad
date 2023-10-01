#ifndef Kernel_Iter_H
#define Kernel_Iter_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * @brief iterative kernel wrapper/(interface) for other kernels
 */
class Kernel_Iter final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter"; }

   private:
    int* istep_ptr;
    int* nstep_ptr;

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iter_H

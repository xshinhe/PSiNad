#ifndef Kernel_Load_DataSet_H
#define Kernel_Load_DataSet_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements a process for loading data in state stucture
 */
class Kernel_Load_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Load_DataSet"; }

   private:
    std::string fn;  ///< filename for loading

    virtual void read_param_impl(Param* PM);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Load_DataSet_H

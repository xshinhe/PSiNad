#ifndef Kernel_Dump_DataSet_H
#define Kernel_Dump_DataSet_H

#include "../core/Kernel.h"

namespace kids {


/**
 * @brief Kernel_Dump_DataSet dump current data with stuctured format
 * @It should be added to the last end (as well as middle) of every solver's
 *  builder
 */
class Kernel_Dump_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Dump_DataSet"; }

   private:
    std::string fn;  ///< filename

    virtual void read_param_impl(Param* PM);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids


#endif  // Kernel_Dump_DataSet_H

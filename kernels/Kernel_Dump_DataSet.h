#ifndef Kernel_Dump_DataSet_H
#define Kernel_Dump_DataSet_H

#include "../core/Kernel.h"

namespace PROJECT_NS {


/**
 * this class implements a process for dumping current data with stuctured format
 * @note
 * this kernel should be appended in the last step of a kernel calling tree (i.e., kernel builder).
 */
class Kernel_Dump_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Dump_DataSet"; }

   private:
    std::string directory;  ///< path for dumping
    std::string fn;         ///< filename (stamp) for dumping
    std::string hdlr_str;   ///< handler type specifier

    virtual void read_param_impl(Param* PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Dump_DataSet_H

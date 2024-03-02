#ifndef Kernel_Load_DataSet_H
#define Kernel_Load_DataSet_H

#include "../core/Kernel.h"

namespace kids {

/**
 * @brief Kernel_Load_DataSet load previous data in state stucture if restart
 *  option is specified
 *  It should be added to the first beginning of every solver's builder.
 */
class Kernel_Load_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Load_DataSet"; }

   private:
    std::string fn;  ///< filename

    virtual void read_param_impl(Param* PM);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids


#endif  // Kernel_Load_DataSet_H

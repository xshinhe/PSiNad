#ifndef Kernel_Load_DataSet_H
#define Kernel_Load_DataSet_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * this class implements a process for loading data in state stucture
 */
class Kernel_Load_DataSet : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    std::string fn;  ///< filename for loading

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Load_DataSet_H

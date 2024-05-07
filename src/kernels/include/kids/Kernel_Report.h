#ifndef Kernel_Report_H
#define Kernel_Report_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Report final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_Report(std::shared_ptr<Kernel_Record> ker) : _recd_ker{ker} { _rules = &(ker->Rules); }

    virtual ~Kernel_Report();

    virtual Status& executeKernel_impl(Status& stat) {
        // clear correlation information
        return 0;
    };

   private:
    std::shared_ptr<Kernel_Record> _recd_ker;
    std::vector<Record_Rule>       _rules;
};

};  // namespace PROJECT_NS


#endif  // Kernel_Report_H

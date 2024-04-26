#ifndef Kernel_Custom_H
#define Kernel_Custom_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

/**
 * This class implements a custom kernel.
 */
class Kernel_Custom final : public Kernel {
   public:
    Kernel_Custom(const std::string& customized_name);

    virtual const std::string getName();

    virtual int getType() const;
};

};  // namespace PROJECT_NS

#endif  // Kernel_Custom_H
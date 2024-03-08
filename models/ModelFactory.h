#ifndef ModelFactory_H
#define ModelFactory_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

extern std::shared_ptr<Kernel> ModelFactory(const std::string& name);

};  // namespace PROJECT_NS


#endif  // ModelFactory_H
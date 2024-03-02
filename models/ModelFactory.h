#ifndef ModelFactory_H
#define ModelFactory_H

#include "../core/Kernel.h"

namespace kids {

extern std::shared_ptr<Kernel> ModelFactory(const std::string& name);

};  // namespace kids


#endif  // ModelFactory_H
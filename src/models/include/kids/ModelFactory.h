#ifndef ModelFactory_H
#define ModelFactory_H

#include "kids/Model.h"

namespace PROJECT_NS {

extern std::shared_ptr<Model> defaultModelFactory(const std::string& name);

};  // namespace PROJECT_NS


#endif  // ModelFactory_H
#ifndef SolverFactory_H
#define SolverFactory_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

extern std::shared_ptr<Kernel> SolverFactory(const std::string& name, std::shared_ptr<Kernel> kmodel);

};  // namespace PROJECT_NS

#endif  // SolverFactory_H

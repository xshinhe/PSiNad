#ifndef SolverFactory_H
#define SolverFactory_H

#include "../core/Kernel.h"

namespace kids {

extern std::shared_ptr<Kernel> SolverFactory(const std::string& name, std::shared_ptr<Kernel> kmodel);

};  // namespace kids

#endif  // SolverFactory_H

#ifndef SolverFactory_H
#define SolverFactory_H

#include "kids/Model.h"
#include "kids/Solver.h"
#include "kids/System.h"

namespace PROJECT_NS {

extern std::shared_ptr<Solver> defaultSolverFactory(const std::string& name, std::shared_ptr<Model> kmodel);

extern std::shared_ptr<Solver> defaultSolverFactory(const std::string& name, std::shared_ptr<System> sys);

};  // namespace PROJECT_NS

#endif  // SolverFactory_H

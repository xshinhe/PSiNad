#include <string>

#include "../utils/definitions.h"
#include "solver.h"

Solver* init_solver(const std::string& solver_name, const Param& parm, Model* pM);

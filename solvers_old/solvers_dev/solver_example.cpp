#include "solver_example.h"

// using namespace ARRAY_EG;

Example_Solver::Example_Solver(Param iparm, Model* pM) : Solver(iparm, pM){};

Example_Solver::~Example_Solver(){};

int Example_Solver::run_impl() { return 0; }

int Example_Solver::run_parallel() { return 0; }

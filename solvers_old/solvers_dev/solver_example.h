#ifndef Example_SOLVER_H
#define Example_SOLVER_H

#include "../solver.h"

class Example_Solver : public Solver {
   public:
    Example_Solver(Param iparm, Model* pM);

    virtual ~Example_Solver();

    static inline std::string name() { return "example"; }

    /**
     * @brief develop you code for the solver
     */
    virtual int run_impl();

    /**
     * @brief develop you code for the solver with mpi parallel
     */
    virtual int run_parallel();

   protected:
    int mydata;
};


#endif  // Example_SOLVER_H

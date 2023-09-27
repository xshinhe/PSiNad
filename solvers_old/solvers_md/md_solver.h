#ifndef MDTraj_Solver_H
#define MDTraj_Solver_H

#include "../solver.h"
#include "traj.h"

class MDTraj_Solver : public Traj_Solver {
   public:
    MDTraj_Solver(Param iparm, Model* pM);
    MDTraj_Solver(const std::string& iparm_str, Model* pM) : MDTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~MDTraj_Solver();

    virtual int ff_calc1(const int& level = 1);

    virtual int traj(TCFnucl& tcfer, const int& N);

    virtual int traj_BAOAB(TCFnucl& tcfer, const int& N);

    virtual int traj_property(const double& dt);

    virtual int sampler(const int& isamp, TCFnucl& tcfer);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

    virtual int run_parallel();

    virtual int init(const int& itraj);

   protected:
    int P, PP, N0, N0N0;
    int pimd_type;
    double bf2;
    std::string pimd_flag;
};


#endif  // MDTraj_Solver_H

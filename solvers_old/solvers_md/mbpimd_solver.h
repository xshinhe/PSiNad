#ifndef MBPIMDTraj_Solver_H
#define MBPIMDTraj_Solver_H

#include "pimd_solver.h"

/**
 * MBPIMDTraj_solver PIMD simulation for fermion, boson and anyon.
 * ref: 10.1103/PhysRevE.106.025309
 *
 */
class MBPIMDTraj_Solver : public PIMDTraj_Solver {
   public:
    MBPIMDTraj_Solver(Param iparm, Model *pM);
    MBPIMDTraj_Solver(const std::string &iparm_str, Model *pM) : MBPIMDTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~MBPIMDTraj_Solver();

    static inline std::string name() { return "mbpimd"; }

    virtual int spring_force();

    virtual int update_p_harm(const num_real &dt);

    virtual int estimator(const int &isamp, TCFnucl &tcfer);

   protected:
    DEFINE_POINTER_PROTECTED(num_real, fV);
    DEFINE_POINTER_PROTECTED(num_real, fE);
    DEFINE_POINTER_PROTECTED(num_real, VHO);
    DEFINE_POINTER_PROTECTED(num_real, dV_spring);
    DEFINE_POINTER_PROTECTED(num_real, dE_spring);
};

#endif  // MBPIMDTraj_Solver_H

#ifndef CMM_SOLVER_H
#define CMM_SOLVER_H
#include "nadtraj.h"

/*
    This shows an example child of NadTraj_Solver: CMM
    All NadTraj solver show detail 3 parts:
    1)  init() function
        initial sampling for mapping variables
    2)  traj() function
        default inherit from (gamma-modified) mean-field traj, (or some revisions)
    3)  correlation(tcf) & sampler() function
        correlation function
*/
class CMM_Solver : public NadTraj_Solver {
   public:
    CMM_Solver(Param iparm, Model* pM);
    CMM_Solver(const std::string& iparm_str, Model* pM) : CMM_Solver(Param::parse(iparm_str), pM){};

    virtual ~CMM_Solver();

    static inline std::string name() { return "cmm"; }

    virtual int init_occ2eac(const int& itraj);

    virtual int init(const int& itraj);

    virtual int kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F);
    virtual int kernel0(num_complex* rhox, const int& F);
    virtual int kernelt(num_complex* rhox, const int& F);

   protected:
    num_real xi0, xit, gamma0, gammat, totact;
};


#endif  // CMM_SOLVER_H

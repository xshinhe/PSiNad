#ifndef WMM_SOLVER_H
#define WMM_SOLVER_H
#include "nadtraj.h"

/*
    This shows an example child of NadTraj_Solver: WMM
    All NadTraj solver show detail 3 parts:
    1)  init() function
        initial sampling for mapping variables
    2)  traj() function
        default inherit from (gamma-modified) mean-field traj, (or some revisions)
    3)  correlation(tcf) & sampler() function
        correlation function
*/
class WMM_Solver : public NadTraj_Solver {
   public:
    WMM_Solver(Param iparm, Model* pM);
    WMM_Solver(const std::string& iparm_str, Model* pM) : WMM_Solver(Param::parse(iparm_str), pM){};

    virtual ~WMM_Solver();

    static inline std::string name() { return "wmm"; }

    virtual int init_occ2eac(const int& itraj);

    virtual int init(const int& itraj);

    virtual int kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F);

    virtual int kernel0(num_complex* rhox, const int& F);
    virtual int kernelt(num_complex* rhox, const int& F);

    virtual int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F);

   protected:
    num_real xi0, xit, gamma0, gammat, totact;
    int adjustp = 0;
};


#endif  // WMM_SOLVER_H

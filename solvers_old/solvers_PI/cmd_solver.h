#ifndef CMD_Solver_H
#define CMD_Solver_H

#include "../solvers_md/pimd_solver.h"

namespace CentroidMDPolicy {
enum _enum { CMD, TRPMD };
const std::map<std::string, _enum> _dict = {{"#CMD", CMD}, {"#RPMD", TRPMD}};
};  // namespace CentroidMDPolicy

class CentroidMD_Solver : public PIMDTraj_Solver {
   public:
    CentroidMD_Solver(Param iparm, Model* pM);
    CentroidMD_Solver(const std::string& iparm_str, Model* pM) : CentroidMD_Solver(Param::parse(iparm_str), pM){};

    virtual ~CentroidMD_Solver();

    static inline std::string name() { return "cmd"; }

    virtual int update_thermo(const num_real& dt);

    //    virtual int run_parallel();

    virtual int init(const int& itraj);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

   protected:
    num_real gam_ad;
    int cmd_type;
    bool constrain_rot;
    std::string CMD_flag;
};
#endif
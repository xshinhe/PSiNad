#ifndef PILDTraj_Solver_H
#define PILDTraj_Solver_H

#include "../../thirdpart/Eigen/Dense"
#include "../../thirdpart/Eigen/SVD"
#include "../solvers_md/pimd_solver.h"

namespace PILDPolicy {
enum _enum { mPILD, PILD };
const std::map<std::string, _enum> _dict = {{"#mPILD", mPILD}, {"#PILD", PILD}};
};  // namespace PILDPolicy

class PILD_Solver : public PIMDTraj_Solver {
   public:
    PILD_Solver(Param iparm, Model* pM);
    PILD_Solver(const std::string& iparm_str, Model* pM) : PILD_Solver(Param::parse(iparm_str), pM){};

    virtual ~PILD_Solver();

    static inline std::string name() { return "pild"; }

    virtual int update_thermo(const num_real& dt);

    virtual int update_p(const num_real& dt);

    //    virtual int run_parallel();

    virtual int init(const int& itraj);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

   protected:
    num_real gam_ad;
    int pild_type;
    bool constrain_rot;
    std::string pild_flag;

    DEFINE_POINTER_PROTECTED(num_real, M_therm);
};

using EigX = Eigen::Matrix<num_real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

int Kabsch(int Natom, int D, num_real* X1, num_real* X2, num_real* m, num_real* R, num_real* r);

#endif
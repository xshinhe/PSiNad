#ifndef MSTRAJ_SOLVER_H
#define MSTRAJ_SOLVER_H

#include "../../models/nad_forcefield/manysite_models.h"
#include "nadtraj.h"

namespace ProductMPSPolicy {
enum _enum { TWA, MMF, CMM };
const std::map<std::string, _enum> _dict = {{"#twa", TWA}, {"#mmf", MMF}, {"#cmm", CMM}};
};  // namespace ProductMPSPolicy

class ProductMPS_Solver : public NadTraj_Solver {
   public:
    ProductMPS_Solver(Param iparm, Model* pM);
    ProductMPS_Solver(const std::string& iparm_str, Model* pM) : ProductMPS_Solver(Param::parse(iparm_str), pM){};

    virtual ~ProductMPS_Solver();

    static inline std::string name() { return "prodmps"; }

    virtual int init(const int& itraj);

    virtual int traj(NAD_TCFer& tcfer, const int& N, const int& F);

    virtual int correlation(const int& isamp, NAD_TCFer& tcfer);

    int M, F, FF, MF, MFF;
    int type;
    num_real scale;

    DEFINE_POINTER(num_complex, rhos);
    DEFINE_POINTER(num_complex, rhos0);
    DEFINE_POINTER(num_complex, rhost);
    DEFINE_POINTER(num_complex, Hs);

    ManySite_ForceField* pForceField;
};


#endif  // MSTRAJ_SOLVER_H
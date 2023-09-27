#ifndef SMM_SOLVER_H
#define SMM_SOLVER_H
#include "nadtraj.h"


class SMM_Solver : public NadTraj_Solver {
   public:
    SMM_Solver(Param iparm, Model* pM);
    SMM_Solver(const std::string& iparm_str, Model* pM) : SMM_Solver(Param::parse(iparm_str), pM){};

    virtual ~SMM_Solver();

    static inline std::string name() { return "smm"; }

    virtual int init_occ2eac(const int& itraj);

    virtual int init(const int& itraj);

    virtual int kernel_scmm(num_complex* rhox, const int& F);

    virtual int kernel0(num_complex* rhox, const int& F);
    virtual int kernelt(num_complex* rhox, const int& F);

   protected:
    DEFINE_POINTER_PROTECTED(num_real, xi1);
};


#endif  // SMM_SOLVER_H

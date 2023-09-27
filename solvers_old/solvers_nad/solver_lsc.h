#ifndef LSC_SOLVER_H
#define LSC_SOLVER_H
#include "nadtraj.h"

class LSC_Solver : public NadTraj_Solver {
   public:
    LSC_Solver(Param iparm, Model* pM);
    LSC_Solver(const std::string& iparm_str, Model* pM) : LSC_Solver(Param::parse(iparm_str), pM){};

    virtual ~LSC_Solver();

    static inline std::string name() { return "lsc"; }

    virtual int init_occ2eac(const int& itraj);

    virtual int init(const int& itraj);

    virtual int kernel_lsc(num_complex* rhox, num_real& gammac, const int& F);

    virtual int kernel0(num_complex* rhox, const int& F);
    virtual int kernelt(num_complex* rhox, const int& F);

   protected:
    num_real gamma0, gammat, sigma0, sigmat, variance;
};


#endif  // LSC_SOLVER_H

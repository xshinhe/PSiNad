#ifndef PMM_Solver_H
#define PMM_Solver_H

#include "multi_nadtraj.h"

class PMM_Solver : public Multi_NadTraj_Solver {
   public:
    PMM_Solver(const Param& iparm, Model* pM);
    PMM_Solver(const std::string& iparm_str, Model* pM) : PMM_Solver(Param::parse(iparm_str), pM){};

    virtual ~PMM_Solver();
    static inline std::string name() { return "pmm"; }

    virtual int traj(NAD_TCFer& tcfer, const int& N, const int& F);
    virtual int multi_evolve_elec(num_complex* Uevolve);
    virtual int multi_ff_calc2();
    virtual int multi_init(const int& itraj);
    virtual int multi_reinit(const int& itraj);
    int virtual kernel0(num_complex* rhox, const int& F);
    int virtual kernelt(num_complex* rhox, const int& F);

   protected:
    num_real scale;
    int isumt;

    DEFINE_POINTER_PROTECTED(num_complex, drhos);
    DEFINE_POINTER_PROTECTED(num_complex, deltarhos);
    DEFINE_POINTER_PROTECTED(num_complex, rhosumt);
    DEFINE_POINTER_PROTECTED(num_complex, rhoavg);
    DEFINE_POINTER_PROTECTED(num_complex, rho_corr);

    int restep;
};

#endif  // PMM_Solver_H

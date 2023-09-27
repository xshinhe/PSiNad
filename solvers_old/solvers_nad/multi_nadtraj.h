#ifndef Multi_NadTraj_H
#define Multi_NadTraj_H

#include "nadtraj.h"

class Multi_NadTraj_Solver : public NadTraj_Solver {
   public:
    Multi_NadTraj_Solver(const Param& iparm, Model* pM, const int& _M = 0);
    Multi_NadTraj_Solver(const std::string& iparm_str, Model* pM) : Multi_NadTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~Multi_NadTraj_Solver();

    static inline std::string name() { return "multinadtraj"; }

    virtual int ff_calc1(const int& level, const bool& refered = false);

    virtual int multi_ff_calc1(const int& level, const bool& refered = false);

    virtual int ff_calc2();

    virtual int multi_ff_calc2();

    virtual int evolve_elec(num_complex* Uevolve);

    virtual int multi_evolve_elec(num_complex* Uevolve);

    virtual int update_p(const num_real& dt_in);

    virtual int update_r(const num_real& dt_in);

    virtual int update_thermo(const num_real& dt_in);

    virtual int init(const int& itraj);

    virtual int multi_init(const int& itraj);

   protected:
    int M, MN, MF, MNN, MFF, MNFF;
    DEFINE_POINTER_PROTECTED(num_real, nrs);
    DEFINE_POINTER_PROTECTED(num_real, nps);
    DEFINE_POINTER_PROTECTED(num_real, nms);
    DEFINE_POINTER_PROTECTED(num_real, nfs);
    DEFINE_POINTER_PROTECTED(num_real, vpeses);
    DEFINE_POINTER_PROTECTED(num_real, grads);
    DEFINE_POINTER_PROTECTED(num_real, hesses);
    DEFINE_POINTER_PROTECTED(num_real, Vs);
    DEFINE_POINTER_PROTECTED(num_real, dVs);
    DEFINE_POINTER_PROTECTED(num_real, ddVs);
    DEFINE_POINTER_PROTECTED(num_real, Es);
    DEFINE_POINTER_PROTECTED(num_real, Ts);
    DEFINE_POINTER_PROTECTED(num_real, dEs);
    DEFINE_POINTER_PROTECTED(num_real, nacvs);
    DEFINE_POINTER_PROTECTED(num_real, Ls);
    DEFINE_POINTER_PROTECTED(num_complex, Ss);
    DEFINE_POINTER_PROTECTED(num_complex, rhos);
    DEFINE_POINTER_PROTECTED(num_complex, Us);
};

#endif  // Multi_NadTraj_H

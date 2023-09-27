#ifndef PPIMDPARA_SOLVER_H
#define PPIMDPARA_SOLVER_H

#include "pimdpara_solver.h"

namespace PIMDPARAHighorderPolicy {
enum _enum { PPI, FD4T3V };
const std::map<std::string, _enum> _dict = {{"ppi", PPI}, {"fd4t3v", FD4T3V}};
};  // namespace PIMDPARAHighorderPolicy

class PPIMDPARATraj_Solver : public PIMDPARATraj_Solver {
   public:
    PPIMDPARATraj_Solver(Param iparm, Model* pM);

    virtual ~PPIMDPARATraj_Solver(){};

    static inline std::string name() { return "ppimd"; }

    virtual int ff_calc1(const int& level = 1);
    virtual int ff_calc2(num_real* trs, num_real* tps, num_real* tvps,
                         // num_real* tgrads,
                         num_real* tfs, const int& level = 1);

    // virtual int update_r(const num_real& dt);

    // virtual int update_p(const num_real& dt);

    // virtual int update_p_harm(const num_real& dt);

    // virtual int caylay_update_half(
    //     const num_real& dt);  // half step of cayley modification is cay(AΔt)^0.5 rather than cay(0.5AΔt)

    // virtual int BAOAB(int& succ,
    //                   int step);  // BAOAB-num in J. Chem. Phys. 145, 024103 (2016) https://doi.org/10.1063/1.4954990

    // virtual int BCOCB(int& succ, int step);

    // virtual int update_thermo(const num_real& dt);

    // virtual int traj(TCFnucl& tcfer, const int& PN);

    // virtual int traj_Middle(TCFnucl& tcfer, const int& PN);

    // virtual int sampler(const int& isamp, TCFnucl& tcfer);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);
    virtual int ppiestimator(const int& isamp, TCFnucl& tcfer);
    virtual int fdestimator(const int& isamp, TCFnucl& tcfer);

    // virtual int run_batch();

    virtual int init(const int& itraj);

    // virtual int final(const int& itraj);

    virtual int all_X2K();
    virtual int all_K2X();
    virtual int all_FX2FK();

    virtual int cal_cent();
    virtual int cal_FDhessian();

   protected:
    int n_term;
    int pon;
    int highorder_type;
    num_real kappa, fppi;
    num_real split_t1, split_v1, split_c1;
    num_real sum_V, sum_Kprim, sum_Kvir, sum_F2, sum_F2v, sum_F2kp, sum_F2kv;
    num_real esti_V, esti_Kprim, esti_Kvir, esti_F2, esti_Bf2;

    int use_liquidne;

    DEFINE_POINTER_PROTECTED(num_real, split_t);
    DEFINE_POINTER_PROTECTED(num_real, split_v);
    DEFINE_POINTER_PROTECTED(num_real, split_c);
    DEFINE_POINTER_PROTECTED(num_real, stag_a);
    DEFINE_POINTER_PROTECTED(num_real, stag_b);
    DEFINE_POINTER_PROTECTED(num_real, stag_c);
    DEFINE_POINTER_PROTECTED(num_real, nhs);
    DEFINE_POINTER_PROTECTED(num_real, nrcent);
};

#endif

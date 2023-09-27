#ifndef SSE_SOLVER_H
#define SSE_SOLVER_H
#include "../../models/nad_forcefield/systembath.h"
#include "../solver.h"
#include "../solvers_nad/nadtcfer.h"

class SSE_Solver : public Solver {
   public:
    enum { SSE_Markovian, SSE_NonMarkovian };
    enum { Tcftype_pop, Tcftype_rho, Tcftype_all };

    const std::map<std::string, int> Map_Tcftype = {{"pop", Tcftype_pop}, {"rho", Tcftype_rho}, {"all", Tcftype_all}};

    SSE_Solver(Param iparm, Model* pM);
    SSE_Solver(const std::string& iparm_str, Model* pM) : SSE_Solver(Param::parse(iparm_str), pM){};

    virtual ~SSE_Solver();

    static inline std::string name() { return "sse"; }

    virtual int init(const int& itraj);

    virtual int correlation(const int& isamp, NAD_TCFer& tcfer);

    virtual int sampler(const int& isamp, NAD_TCFer& tcfer);

    virtual int action_on_wavafunction(num_complex* eacnew, num_complex* eac, const num_real& t);

    virtual int action_on_wavafunction_adia(num_complex* eacnew, num_complex* eac, const num_real& t);

    virtual int time_kernel(num_complex* Ker, const num_real& t);

    virtual int call_memory_array_adia(const num_real& dt);

    virtual int sse(NAD_TCFer& tcfer, const int& N, const int& F);  // gamma-modified mean-field trajectory

    virtual int sse_adia(NAD_TCFer& tcfer, const int& N, const int& F);  // gamma-modified mean-field trajectory

    virtual int action_on_wavafunction_adia_mem(num_complex* eacnew, num_complex* eacold, const int& hstep,
                                                const bool& hupdate);

    virtual int run_parallel();

   protected:
    int N, F, NN, FF, NF, NFF, NNF, NNFF, FFFF, L;
    int nbath, Nb, NbFF;
    int istep, sstep, nstep, isamp, nsamp, itraj, ntraj, ispec, nspec;
    int occ0, occt;
    num_real vscale;
    num_real dt, tend, beta;

    DEFINE_POINTER_PROTECTED(num_real, rand1);
    DEFINE_POINTER_PROTECTED(num_real, rand2);
    DEFINE_POINTER_PROTECTED(num_real, Eele);
    DEFINE_POINTER_PROTECTED(num_real, DE);
    DEFINE_POINTER_PROTECTED(num_real, w);
    DEFINE_POINTER_PROTECTED(num_real, tanhqwb);
    DEFINE_POINTER_PROTECTED(num_real, xicoeff1);
    DEFINE_POINTER_PROTECTED(num_real, xicoeff2);
    DEFINE_POINTER_PROTECTED(num_real, alpha_pref);
    DEFINE_POINTER_PROTECTED(num_real, CL);
    DEFINE_POINTER_PROTECTED(num_real, DEpW);
    DEFINE_POINTER_PROTECTED(num_real, DEmW);
    DEFINE_POINTER_PROTECTED(num_real, coeff_DEpW);
    DEFINE_POINTER_PROTECTED(num_real, coeff_DEmW);
    DEFINE_POINTER_PROTECTED(num_real, Hele);
    DEFINE_POINTER_PROTECTED(num_real, Tele);
    DEFINE_POINTER_PROTECTED(num_real, T0);
    DEFINE_POINTER_PROTECTED(num_real, Xzero);
    DEFINE_POINTER_PROTECTED(num_real, Xtran);
    DEFINE_POINTER_PROTECTED(num_real, Qzero);
    DEFINE_POINTER_PROTECTED(num_real, Qtran);

    DEFINE_POINTER_PROTECTED(num_complex, BL);
    DEFINE_POINTER_PROTECTED(num_complex, tmpbarr1);
    DEFINE_POINTER_PROTECTED(num_complex, tmpbarr2);

    DEFINE_POINTER_PROTECTED(num_complex, eac);
    DEFINE_POINTER_PROTECTED(num_complex, eac_adia);
    DEFINE_POINTER_PROTECTED(num_complex, eac0);
    DEFINE_POINTER_PROTECTED(num_complex, rho0);
    DEFINE_POINTER_PROTECTED(num_complex, rhot);
    DEFINE_POINTER_PROTECTED(num_complex, fact1);
    DEFINE_POINTER_PROTECTED(num_complex, fact2);
    DEFINE_POINTER_PROTECTED(num_complex, fact3);
    DEFINE_POINTER_PROTECTED(num_complex, fact4);
    DEFINE_POINTER_PROTECTED(num_complex, eactmp);
    DEFINE_POINTER_PROTECTED(num_complex, eactmp1);
    DEFINE_POINTER_PROTECTED(num_complex, eactmp2);
    DEFINE_POINTER_PROTECTED(num_complex, eacadd1);
    DEFINE_POINTER_PROTECTED(num_complex, eacadd2);
    DEFINE_POINTER_PROTECTED(num_complex, eactran);
    DEFINE_POINTER_PROTECTED(num_complex, crand1);
    DEFINE_POINTER_PROTECTED(num_complex, crand2);
    DEFINE_POINTER_PROTECTED(num_complex, invexpiEt);
    DEFINE_POINTER_PROTECTED(num_complex, U);
    DEFINE_POINTER_PROTECTED(num_complex, expiwt);
    DEFINE_POINTER_PROTECTED(num_complex, expiDEt);
    DEFINE_POINTER_PROTECTED(num_complex, invexpiwt);
    DEFINE_POINTER_PROTECTED(num_complex, invexpiDEt);
    DEFINE_POINTER_PROTECTED(num_complex, expiwdht);
    DEFINE_POINTER_PROTECTED(num_complex, invexpiwdht);
    DEFINE_POINTER_PROTECTED(num_complex, expiwt_now);
    DEFINE_POINTER_PROTECTED(num_complex, invexpiwt_now);
    DEFINE_POINTER_PROTECTED(num_complex, Xtmpj);
    DEFINE_POINTER_PROTECTED(num_complex, Xtmp2j);
    DEFINE_POINTER_PROTECTED(num_complex, xi);
    DEFINE_POINTER_PROTECTED(num_complex, time_ker);
    DEFINE_POINTER_PROTECTED(num_complex, mem_arr);

    SystemBath_ForceField* pForceField;

    int tcf_type;
    int ini_type;
    std::string tcf_flag;
    std::string ini_flag;
};

#endif  // SSE_SOLVER_H

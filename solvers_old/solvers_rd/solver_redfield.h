#ifndef REDFIELD_SOLVER_H
#define REDFIELD_SOLVER_H

#include "../../models/nad_forcefield/systembath.h"
#include "../solver.h"

class RedField_Solver : public Solver {
   public:
    RedField_Solver(Param iparm, Model* pM);
    RedField_Solver(const std::string& iparm_str, Model* pM) : RedField_Solver(Param::parse(iparm_str), pM){};

    virtual ~RedField_Solver();

    static inline std::string name() { return "redfield"; }

    virtual int run_impl();

    int action_on_wavafunction(num_complex* rhonew, num_complex* rhoold, const num_real& t);

    int calc_R_tensor();

   protected:
    int N, F, NN, FF, NF, NFF, NNF, NNFF, FFFF, L;
    int nbath, Nb, NbFF;
    int istep, sstep, nstep, isamp, nsamp, itraj, ntraj, ispec, nspec;
    int occ0, occt;
    num_real vscale;
    num_real dt, tend;

    DEFINE_POINTER_PROTECTED(num_real, Eele);
    DEFINE_POINTER_PROTECTED(num_real, Tele);
    DEFINE_POINTER_PROTECTED(num_real, Hele);
    DEFINE_POINTER_PROTECTED(num_real, mDE);
    DEFINE_POINTER_PROTECTED(num_real, Qtran);

    DEFINE_POINTER_PROTECTED(num_complex, C_mDE);
    DEFINE_POINTER_PROTECTED(num_complex, GM_tensor);
    DEFINE_POINTER_PROTECTED(num_complex, R_tensor);

    DEFINE_POINTER_PROTECTED(num_complex, eac0);
    DEFINE_POINTER_PROTECTED(num_complex, rho0);
    DEFINE_POINTER_PROTECTED(num_complex, rhoadia);
    DEFINE_POINTER_PROTECTED(num_complex, rhodia);

    // for RK-4
    DEFINE_POINTER_PROTECTED(num_complex, rho1);
    DEFINE_POINTER_PROTECTED(num_complex, rho2);
    DEFINE_POINTER_PROTECTED(num_complex, rho3);
    DEFINE_POINTER_PROTECTED(num_complex, rho4);
    DEFINE_POINTER_PROTECTED(num_complex, rhotmp);

    SystemBath_ForceField* pForceField;
};

#endif  // REDFIELD_SOLVER_H

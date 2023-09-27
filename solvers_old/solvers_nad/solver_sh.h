#ifndef Hopping_Solver_H
#define Hopping_Solver_H

#include "../../utils/definitions.h"
#include "../../utils/nad_utils.h"
#include "nadtraj.h"

/**
 * realized methods:
 * 0) FSSH                                  | checked
 * 1) PC-FSSH: doi:10.1063/1.3603447        | [init]
 * 3) GFSH: doi:10.1063/1.4935971           | [init]
 * 4) NDM: doi:10.1063/1.1368388            | [init]
 * 4) SCDM: doi:10.1063/1.1648306           | [init]
 * 5) AFSSH: doi:10.1063/1.3506779          | [init]
 * 6) DISH: doi:10.1063/1.4757100           | [init]
 */

namespace hopping {
enum _enum { FSSH, DISH, GFSH, AFSSH, NDM, SCDM };
const std::map<std::string, _enum> _dict = {
    {"fssh", FSSH}, {"dish", DISH}, {"gfsh", GFSH}, {"afssh", AFSSH}, {"ndm", NDM}, {"scdm", SCDM},
};
};  // namespace hopping

class Hopping_Solver : public NadTraj_Solver {
   public:
    Hopping_Solver(Param iparm, Model* pM);
    Hopping_Solver(const std::string& iparm_str, Model* pM) : Hopping_Solver(Param::parse(iparm_str), pM){};

    virtual ~Hopping_Solver();

    static inline std::string name() { return "hopping"; }

    virtual int init_occ2eac(const int& itraj);

    virtual int init(const int& itraj);

    virtual inline num_real nucl_Ekin() {
        num_real res = 0;
        for (int j = 0; j < N; ++j) res += 0.5f * np[j] * np[j] / nm[j];
        return res;
    }

    virtual int SE_Hamiltonian(const bool& refered);

    virtual int phase_correction();

    virtual int hopping_choose(const int& iocc, num_complex* rhox, num_complex* H, num_real& dt);

    virtual int hopping_direction(num_real* direction, const int& to);

    virtual int hopping_impulse(num_real* direction, const int& occ_to);

    virtual int time_calc();

    virtual int coherent_evolve(const num_real& dt);

    virtual int decoherent_evolve(const num_real& dt);

    virtual int hopping_collapse(const num_real& dt, const int& to);

    virtual int ff_calc1(const int& level, const bool& refered = false);

    virtual int kernel_fssh(num_complex* rhox, const int& F);

    virtual int kernel0(num_complex* rhox, const int& F);

    virtual int kernelt(num_complex* rhox, const int& F);

    virtual int check_break(int& succ);

    virtual int ff_calc2();

    virtual int traj(NAD_TCFer& tcfer, const int& N, const int& F);

   protected:
    bool reflect;
    bool pcorrect;
    bool terminate;
    bool dish_prefer_deh1;
    bool dish_prefer_deh2;
    bool dish_project_self;

    DEFINE_POINTER_PROTECTED(num_complex, arg_nr);  ///<?
    DEFINE_POINTER_PROTECTED(num_complex, arg_np);

    num_real dishw;

    DEFINE_POINTER_PROTECTED(num_real, time_coh);
    DEFINE_POINTER_PROTECTED(num_real, time_los);
    DEFINE_POINTER_PROTECTED(num_real, dnp_diss);
    DEFINE_POINTER_PROTECTED(num_real, probs_dish);
    DEFINE_POINTER_PROTECTED(int, idx_dish);
    DEFINE_POINTER_PROTECTED(num_real, direction);
    DEFINE_POINTER_PROTECTED(num_complex, drho_diss);
    DEFINE_POINTER_PROTECTED(num_complex, Uh);

    int Fadd1;

    int hopping_strategy;
};



#endif  // Hopping_Solver_H

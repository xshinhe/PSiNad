#ifndef NADTRAJ_H
#define NADTRAJ_H

#include <fstream>
#include <iostream>
#include <sstream>

#include "../../utils/types.h"
#include "../solvers_md/traj.h"
#include "nadtcfer.h"

class NadTraj_Solver : public Traj_Solver {
   public:
    NadTraj_Solver(const Param& iparm, Model* pM);
    NadTraj_Solver(const std::string& iparm_str, Model* pM) : NadTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~NadTraj_Solver();

    static inline std::string name() { return "nadtraj"; }

    virtual int ff_calc1(const int& level, const bool& refered = false);

    virtual int ff_calc2();

    virtual int evolve_elec(num_complex* Uevolve);

    virtual int init_occ2eac(const int& itraj);

    virtual int init_ofs(const int& itraj);

    virtual int init(const int& itraj);

    virtual int final(const int& itraj);

    virtual int check_break(int& succ);

    virtual int rst_output(const int& itraj);

    virtual int rst_read(const int& itraj);

    virtual int traj_property(const num_real& dt);

    virtual int traj(NAD_TCFer& tcfer, const int& N, const int& F);

    virtual int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F);

    int sampler(const int& isamp, NAD_TCFer& tcfer);

    virtual int kernel0(num_complex* rhox, const int& F);

    virtual int kernelt(num_complex* rhox, const int& F);

    virtual int correlation(const int& isamp, NAD_TCFer& tcfer);

    virtual int run_impl();

    virtual int run_parallel();

   protected:
    int F, FF, NF, NFF, NNF, NNFF, FFFF;

    int occ0, occt;

    DEFINE_POINTER_PROTECTED(num_real, mvc);

    DEFINE_POINTER_PROTECTED(num_real, fmean);  ///< electronic weighted mean force
    DEFINE_POINTER_PROTECTED(num_real, fcorr);  ///< correction to force

    DEFINE_POINTER_PROTECTED(num_real, V);
    DEFINE_POINTER_PROTECTED(num_real, dV);
    DEFINE_POINTER_PROTECTED(num_real, ddV);
    DEFINE_POINTER_PROTECTED(num_real, E);
    DEFINE_POINTER_PROTECTED(num_real, T);
    DEFINE_POINTER_PROTECTED(num_real, T0);
    DEFINE_POINTER_PROTECTED(num_real, dE);
    DEFINE_POINTER_PROTECTED(num_real, ddE);
    DEFINE_POINTER_PROTECTED(num_real, L);
    DEFINE_POINTER_PROTECTED(num_real, nacv);

    DEFINE_POINTER_PROTECTED(num_complex, eac0);
    DEFINE_POINTER_PROTECTED(num_complex, rho0);
    DEFINE_POINTER_PROTECTED(num_complex, rhot);
    DEFINE_POINTER_PROTECTED(num_complex, U);
    DEFINE_POINTER_PROTECTED(num_complex, H);
    DEFINE_POINTER_PROTECTED(num_complex, dH);
    DEFINE_POINTER_PROTECTED(num_complex, ddH);
    DEFINE_POINTER_PROTECTED(num_complex, S);
    DEFINE_POINTER_PROTECTED(num_complex, dL);
    DEFINE_POINTER_PROTECTED(num_complex, ddL);
    DEFINE_POINTER_PROTECTED(num_complex, eac);
    DEFINE_POINTER_PROTECTED(num_complex, eacf);
    DEFINE_POINTER_PROTECTED(num_complex, eacb);
    DEFINE_POINTER_PROTECTED(num_complex, gmat);
    DEFINE_POINTER_PROTECTED(num_complex, rho);

    Nad_ForceField* pForceField;

    num_complex tcf_weight = phys::math::iu;

    int mem_type;
    int rep_type;
    int eom_type;
    int tcf_type;
    int ini_type;
    int dyn_type;

    std::ofstream ofs_ESAMP;
    std::ofstream ofs_ETRAJ;
};



#endif  // NADTRAJ_H

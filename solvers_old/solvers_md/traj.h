#ifndef TRAJ_H
#define TRAJ_H

#include "../../models/forcefieldbase.h"
#include "../../models/thermostat/thermostat.h"
#include "../solver.h"

using TCFnucl = num_real;

/**
 * @brief Traj_Solver: abstract interface for trajectory-based approach
 * NOT use it as MD solver!!!
 */
class Traj_Solver : public Solver {
   public:
    Traj_Solver(const Param& iparm, Model* pM);
    Traj_Solver(const std::string& iparm_str, Model* pM) : Traj_Solver(Param::parse(iparm_str), pM){};

    virtual ~Traj_Solver();

    static inline std::string name() { return "traj"; }

    virtual int ff_calc1(const int& level);

    virtual int init_ofs(const int& itraj);

    virtual int init(const int& itraj);  // by default or from hdf5 file

    virtual int final(const int& itraj);  // write hdf5 restart file and something else

    virtual int rst_output(const int& traj_in);

    virtual int rst_read(const int& traj_in);

    virtual int check_break(int& succ);

    virtual int traj_property(const num_real& dt);

    virtual int update_r(const num_real& dt);

    virtual int update_p(const num_real& dt);

    virtual int update_thermo(const num_real& dt);

    virtual int traj(TCFnucl& tcfer, const int& N);

    virtual int traj_velocityverlet(TCFnucl& tcfer, const int& N);

    virtual int sampler(const int& isamp, TCFnucl& tcfer);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

    virtual int run_impl();

    virtual int run_parallel();

   protected:
    int N, NN;

    int level;
    int Ndim;
    int Natom;

    // for monte carlo
    int ntraj, itraj;
    int nstep, istep, sstep;
    int nsamp, isamp, lsamp;
    int nspec, ispec;
    int nrespa;

    // time-step setting for dynamics
    num_real tstart = 0.0f, tend;  //< @TODO tstart to be used
    num_real dt, halfdt, smalldt, largedt;

    num_real Hhere, Lhere, Shere, Khere, Vhere;

    DEFINE_POINTER_PROTECTED(num_real, Htot);  ///< Hamiltonian time sequence
    DEFINE_POINTER_PROTECTED(num_real, Ltot);  ///< Lagrangiain time sequence
    DEFINE_POINTER_PROTECTED(num_real, Stot);  ///< action time sequence
    DEFINE_POINTER_PROTECTED(num_real, Ktot);  ///< kinetic energy time sequence
    DEFINE_POINTER_PROTECTED(num_real, Vtot);  ///< potential energy time sequence

    DEFINE_POINTER_PROTECTED(num_real, nr0);  ///< initial sampling of nuclear position
    DEFINE_POINTER_PROTECTED(num_real, np0);  ///< initial sampling of nuclear momentum

    DEFINE_POINTER_PROTECTED(num_real, nr);  ///< (dynamical) position DOFs
    DEFINE_POINTER_PROTECTED(num_real, np);  ///< (dynamical) momentum DOFs
    DEFINE_POINTER_PROTECTED(num_real, nm);  ///< mass DOFs
    DEFINE_POINTER_PROTECTED(num_real, nf);  ///< (dynamical) gradient DOFs (negative force)

    DEFINE_POINTER_PROTECTED(num_real, vpes);  ///< saves potential energy from pForceField
    DEFINE_POINTER_PROTECTED(num_real, grad);  ///< saves gradient from pForceField
    DEFINE_POINTER_PROTECTED(num_real, hess);  ///< saves hessain from pForceField

    Thermostat* pThermo;         ///< bind with thermostat
    BO_ForceField* pForceField;  ///< bind with forceField

    std::ofstream ofs_SAMP;  ///< ofstream for initial sampling
    std::ofstream ofs_TRAJ;  ///< ofstream for traj information
    std::ofstream ofs_ENER;  ///< ofstream for energy information
};



#endif  // TRAJ_H

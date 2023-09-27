#ifndef QCPI_SOLVER_H
#define QCPI_SOLVER_H
#include "nadtraj.h"

/*
    This shows an example child of NadTraj_Solver: CMM
    All NadTraj solver show detail 3 parts:
    1)  init() function
        initial sampling for mapping variables
    2)  traj() function
        default inherit from (gamma-modified) mean-field traj, (or some revisions)
    3)  correlation(tcf) & sampler() function
        correlation function
*/
class QCPI_Solver : public NadTraj_Solver {
   public:
    QCPI_Solver(Param iparm, Model* pM);
    QCPI_Solver(const std::string& iparm_str, Model* pM) : QCPI_Solver(Param::parse(iparm_str), pM){};

    virtual ~QCPI_Solver();

    static inline std::string name() { return "qcpi"; }

    inline int get_index(int* arr, const int& num, const int& base, const int& len) {
        plFunction();
        // from lower => higher site
        for (int i = 0, c = num; i < len; ++i) arr[i] = c % base, c /= base;
        return 0;
    }

    inline int get_Nnum(int* arr, const int& base, const int& len) {
        plFunction();
        int Nnum = 0;
        for (int i = 0, wt = 1; i < len; ++i) Nnum += arr[i] * wt, wt *= base;
        return Nnum;
    }

    virtual int propagate_tensor();

    virtual num_complex propagate_nuc(num_real* nr2_trj, num_real* np2_trj,  // updated traj
                                      num_real* nr1_trj, num_real* np1_trj,  // old traj
                                      const int& prev, const int& next, const num_real& dt);
    virtual int init(const int& itraj);

    virtual int kernel0(num_complex* rhox, const int& F);

    virtual int kernelt(num_complex* rhox, const int& F);

    virtual int traj(NAD_TCFer& tcfer, const int& N, const int& F);

   protected:
    int Lx, Ly, Nx, Ny, mem, niter, Lx_now, Ly_now;

    DEFINE_POINTER_PROTECTED(int, idx_arr);
    DEFINE_POINTER_PROTECTED(num_real, nrs);
    DEFINE_POINTER_PROTECTED(num_real, nps);
    DEFINE_POINTER_PROTECTED(num_real, nrs_copy);
    DEFINE_POINTER_PROTECTED(num_real, nps_copy);
    DEFINE_POINTER_PROTECTED(num_complex, nx);
    DEFINE_POINTER_PROTECTED(num_complex, ny);
    DEFINE_POINTER_PROTECTED(num_complex, nx_copy);
    DEFINE_POINTER_PROTECTED(num_complex, ny_copy);
};


#endif  // QCPI_SOLVER_H

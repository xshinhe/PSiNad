#ifndef HEOM_SOLVER_H
#define HEOM_SOLVER_H

#include <omp.h>

#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../models/nad_forcefield/systembath.h"
#include "../../utils/definitions.h"
#include "../../utils/sparse_utils.h"
#include "../solver.h"

class HEOM_Solver : public Solver {
   public:
    enum { BasisExpan_Matsubara, BasisExpan_Fitting, BasisExpan_Closure };
    enum { TCF_rtype, TCF_itype, TCF_xtype };

    HEOM_Solver(Param iparm, Model* pM);

    HEOM_Solver(const std::string& iparm_str, Model* pM) : HEOM_Solver(Param::parse(iparm_str), pM){};

    virtual ~HEOM_Solver();

    static inline std::string name() { return "heom"; }

    virtual int Init_Bath(const int& ib, const int& Nr, const int& Ni, num_real& gdph_m, num_complex* tcf_coef,
                          num_complex* tcf_zero, num_complex* tcf_deri);

    virtual int Generate_ADOs();

    virtual int run_impl();

    virtual int run_parallel() { return run_impl(); };

    virtual int Evolve(num_complex* sigma_tot_t1, num_complex* sigma_tot_t2);

    virtual bool if_filteration(int* arr);

    long long int nchoosek(unsigned int n, unsigned int k);

    long long int findpos(int* arr, const int& N);

    int findarr(int* arr, const long long int& iX, const int& N);

    int N, H, F, FF, nbath;
    int Nexpan_M, Nexpan_Nr, Nexpan_Ni;
    int Nchs;
    long long int** CNK;

    int NADOs, Nconn;

    DEFINE_POINTER(int, csr_iADOs);
    DEFINE_POINTER(int, csr_connec_LD);
    DEFINE_POINTER(int, csr_connec);
    DEFINE_POINTER(int, csr_type);
    DEFINE_POINTER(int, csr_ivalue);
    DEFINE_POINTER(num_complex, csr_cvalue);

    DEFINE_POINTER(int, Tcftype);
    DEFINE_POINTER(int, tcf_site);

    DEFINE_POINTER(num_real, Esys);
    DEFINE_POINTER(num_complex, tcf_coef);
    DEFINE_POINTER(num_complex, tcf_zero);
    DEFINE_POINTER(num_complex, tcf_deri);
    DEFINE_POINTER(num_complex, sigma_tot);
    DEFINE_POINTER(num_complex, rdm_arr);
    DEFINE_POINTER(num_complex, Tsys);
    DEFINE_POINTER(num_complex, Hsys);
    DEFINE_POINTER(num_complex, T0);
    DEFINE_POINTER(num_complex, Hsys_time);
    DEFINE_POINTER(num_complex, eac0);
    DEFINE_POINTER(num_complex, rho0);

    SparseMat<num_complex>**LQ_list, **SQ_list, **LQLQ_list, *ALQLQ, *LHsys, *LHsysDphs;

    int BasisExpanType;
    std::string basis_file;
    SystemBath_ForceField* pForceField;

    // for running
    int sstep, nstep, adia_pure;
    int occ0;
    num_real dt, tend, vscale;
};

#endif  // HEOM_SOLVER_H

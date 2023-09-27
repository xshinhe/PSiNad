#ifndef SCF_Solver_H
#define SCF_Solver_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../solver.h"
#include "basis_set.h"

struct Atom {
    int znum;
    double x, y, z;
};

class SCF_Solver : public Solver {
   public:
    SCF_Solver(const Param &iparm, Model *pM);
    SCF_Solver(const std::string &iparm_str, Model *pM) : SCF_Solver(Param::parse(iparm_str), pM){};

    virtual ~SCF_Solver();

    static inline std::string name() { return "scf"; }

    int print_info();
    int call_ao_integral();
    int call_ao_integral_simple();

#ifdef USE_LIBCINT
    int call_ao_integral_libcint();
#endif

    int initial_guess();
    int solve_orthonormalize();
    int iter(const int &cnt);
    int iter_simple(const int &cnt);
    int iter_diis(const int &cnt);
    int build_Fock();
    int solve_Fock();
    int sync_mo2ao();
    int analysis(const int &cnt);
    int scf_report(const int &cnt);

    virtual int run_impl();

    virtual int run_parallel() { return run_impl(); };

    double E_nucl_repul();


    Atom *atoms;  ///< save atoms

    std::vector<Atomic_BasisSet> basis_list;  ///< save Atomic_Basis for each kind of atoms
    // Atomic_BasisSet** basis_pp_unique_list;
    Atomic_BasisSet **ppbasis_helper;

    std::map<int, int> atom_count;   ///< count No. of each kind of atoms
    std::map<int, int> atom_bslink;  ///< mapping znum to Atomic_Basis (index in basis_list)

    double *S_AO;    ///< overlap intergal in AOs
    double *H_AO;    ///< single-electron intergal in AOs
    double *ERI_AO;  ///< double-electrons (ERI) intergal in AOs

    double *matX;  ///< S^(-1/2)
    double *matY;  ///< S^(+1/2)

    double *G_AO;  ///< double-electron maritx (functional of MOs) in AOs
    double *J_AO;  ///< comloub intergal maritx (functional of MOs) in AOs
    double *K_AO;  ///< exchange integral maritx (functional of MOs) in AOs

    double *F_AO;       ///< fock matrix in AOs
    double *C_AO;       ///< coeffients matrix in AOs
    double *P_AO;       ///< population matrix in AOs
    double *P_AO_last;  ///< population matrix in AOs in last step

    double *F_MO;  ///< fock matrix in MOs (only orthonormalized but not canonical MOs)
    double *C_MO;  ///< coeffients matrix in MOs
    double *P_MO;  ///< population matrix in MOs
    double *E_MO;  ///< eigenvalues in MOs

    double *P_NO;  ///< population matrix in canonical (natural) MOs

    int maxiter;
    double prec;

    int Q;      ///< Charge
    int S;      ///< Spin multiplicity
    int Nelec;  ///< No. of electrons
    int Nocc;   ///< No. of occupations
    int Natom;  ///< No. of atoms
    int Nb;     ///< No. of basis
    int NNb;
    int NNNb;

    double E, E_last;
};


#endif  // SCF_Solver_H

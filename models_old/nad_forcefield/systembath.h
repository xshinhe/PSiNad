#ifndef SystemBath_H
#define SystemBath_H

#include "../bath/bath.h"
// #include "../forcefieldbase.h"

namespace SystemBathPolicy {
enum _enum_Hsys {
    SpinBosonHsys,
    FMOHsys,
    SF3aHsys,
    SF3bHsys,
    SF3cHsys,
    SF5aHsys,
    SF5bHsys,
    FCPHsys,
    AGGHsys,
    CYCHsys,
    GivenHsys,
};
const std::map<std::string, int> _dict_Hsys = {
    {"#sb", SpinBosonHsys},  // sigma_z is coulped
    {"#fmo", FMOHsys},       //
    {"#sf3a", SF3aHsys},     //
    {"#sf3b", SF3bHsys},     //
    {"#sf5a", SF5aHsys},     //
    {"#sf5b", SF5bHsys},     //
    {"#fcp", FCPHsys},       //
    {"#agg", AGGHsys},       //
    {"#cyc", CYCHsys},       //
    {"#rub", CYCHsys},       //
    {"read", GivenHsys}      // read from "Hsys.dat"
};
};  // namespace SystemBathPolicy

class SystemBath_ForceField final : public Model {
   public:
    static inline std::string name() { return "systembath"; }

    static inline double J_Ohmic(const double& w);

    static inline double J_Debye(const double& w);

    static inline double J_Refit(const double& w, double* W_arr, double* C_arr, const int& Nb);

    static inline double J(const double& w, double* W_arr = nullptr, double* C_arr = nullptr, const int& Nb = 0);

    static inline double fun_Cw_Re(const double& w, double* W_arr, double* C_arr, const int& Nb);

    static inline int fun_Cw(num_complex* Cw_arr, double* w, const int& Nw, double* W_arr, double* C_arr,
                             const int& Nb);

   private:
    // parameters
    int N;
    int nbath;
    int Nb;
    int F;
    int NN, NF, FF, NFF, NNFF;
    int L;  // no. of nonzero variables in each Q

    // references
    num_real* Hsys;    /* Hamiltonian for system part */
    num_real* Q;       /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */
    num_real* omegas;  ///< save discrete frequencies (only for simple model, L=1)
    num_real* coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    num_real* CL;      ///< save coupling coefficients with Qj (Qj has L no. of nonzero elements)
    num_real* QL;      ///< save coulping matrix, each and L no. of nonzero elements
    num_real* Xnj;     ///< used in Stochastic Schrodinger Equation Methods

    // options
    Option Coup_option;
    Option Hsys_option;
    Option Lamb_option;
    Option Spec_option;


    virtual void read_param_impl(Param* P);
    virtual void init_data_impl(State* S);
    virtual int exec_kernel_impl(int stat = -1);

    DEFINE_POINTER(num_real, Hsys); /* Hamiltonian for system part */
    DEFINE_POINTER(num_real, Q);    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */

    Bath* mybath; /* description for bath */

    DEFINE_POINTER(num_real, omegas);  ///< save discrete frequencies (only for simple model, L=1)
    DEFINE_POINTER(num_real, coeffs);  ///< save coupling coefficients (only for simple model, L=1)

    // for sparse representation of CQ
    int L;                          // no. of nonzero variables in each Q
    DEFINE_POINTER(num_real, CL);   ///< save coupling coefficients with Qj (Qj has L no. of nonzero elements)
    DEFINE_POINTER(num_real, QL);   ///< save coulping matrix, each and L no. of nonzero elements
    DEFINE_POINTER(num_real, Xnj);  ///< (used in Stochastic Schrodinger Equation Methods)

   protected:
    int Nb, nbath, coup_type;
    bool first_call = true;

    int Hsys_type;
    std::string Hsys_flag;
};

#endif  // SystemBath_H

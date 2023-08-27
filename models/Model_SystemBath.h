#ifndef SystemBath_H
#define SystemBath_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Model_Bath.h"

namespace PROJECT_NS {

DEFINE_POLICY(SystemPolicy,  //
              SB,            //
              FMO,           //
              SF3a,          //
              SF3b,          //
              SF3c,          //
              SF5a,          //
              SF5b,          //
              FCP,           //
              AGG,           //
              CYC,           //
              Read);         //

DEFINE_POLICY(CouplingPolicy,  //
              SB,              //
              SE,              //
              Read);           //

DEFINE_POLICY(NSampPolicy,
              Wigner,     //
              Classical,  //
              QCT);

class Model_SystemBath final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_SystemBath"; }

    Model_SystemBath() {
        push(std::shared_ptr<Model_Bath>(new Model_Bath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of bath
    int Nb;     // discrete no.
    int L;      // no. of nonzero variables in each Q

    // system & coupling
    num_real* Hsys; /* Hamiltonian for system part */
    num_real* Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */
    num_real* CL;   ///< save coupling coefficients with Qj (Qj has L no. of nonzero elements)
    num_real* QL;   ///< save coulping matrix, each and L no. of nonzero elements
    num_real* Xnj;  ///< used in Stochastic Schrodinger Equation Methods

    // bath
    num_real* omegas;  ///< save discrete frequencies (only for simple model, L=1)
    num_real* coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    num_real* x_sigma;
    num_real* p_sigma;

    // integrator
    num_real *x, *p, *m;

    // model
    num_real* mass;
    num_real *vpes, *grad, *hess;
    num_real *V, *dV, *ddV;

    // options
    SystemPolicy::_type system_type;
    BathPolicy::_type bath_type;
    CouplingPolicy::_type coupling_type;
    NSampPolicy::_type nsamp_type;

    virtual void read_param_impl(Param* P);
    virtual void init_data_impl(DataSet* DS);
    virtual void init_calc_impl(int stat = -1);
    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // SystemBath_H

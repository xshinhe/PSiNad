#ifndef SystemBath_H
#define SystemBath_H

#include "kids/Kernel.h"
#include "kids/Model_Bath.h"
#include "kids/Policy.h"

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
    virtual const std::string getName();

    virtual int getType() const;

    Model_SystemBath() {
        appendChild(std::shared_ptr<Model_Bath>(new Model_Bath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of bath
    int Nb;     // discrete no.
    int L;      // no. of nonzero variables in each Q

    // system & coupling
    kids_real* Hsys; /* Hamiltonian for system part */
    kids_real* Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */
    kids_real* CL;   ///< save coupling coefficients with Qj (Qj has L no. of nonzero elements)
    kids_real* QL;   ///< save coulping matrix, each and L no. of nonzero elements
    kids_real* Xnj;  ///< used in Stochastic Schrodinger Equation Methods

    // bath
    kids_real* omegas;  ///< save discrete frequencies (only for simple model, L=1)
    kids_real* coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    kids_real* x_sigma;
    kids_real* p_sigma;

    // integrator
    kids_real *x, *p, *m;

    // model
    kids_real* mass;
    kids_real *vpes, *grad, *hess;
    kids_real *V, *dV, *ddV;

    // options
    SystemPolicy::_type   system_type;
    BathPolicy::_type     bath_type;
    CouplingPolicy::_type coupling_type;
    NSampPolicy::_type    nsamp_type;

    virtual void    setInputParam_impl(std::shared_ptr<Param>& PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet>& DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // SystemBath_H

#ifndef SystemBath_H
#define SystemBath_H

#include "kids/Model.h"
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

class Model_SystemBath final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_SystemBath() {
        appendChild(std::shared_ptr<Model_Bath>(new Model_Bath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of baths
    int Nb;     // no. of discrete modes for a bath
    int L;      // no. of nonzero variables in each interaction matrix Q

    // system & coupling
    kids_real* Hsys;  ///< [F * F] (electonic) system Hamiltonian matrix
    kids_real* Kmat;  ///< [N * N] (nuclear) oscillation strength matrix
    kids_real* Qmat;  ///< [N * F * F] coupling matrix
    kids_real* Q;     ///< [nbath * FF] interaction matrix for different bath
    kids_real* CL;    ///< [L * Nb] discretized coefficients times nonzero terms of Qj
    kids_real* QL;    ///< [L * Nb * F * F] save coulping matrix, each and L no. of nonzero elements
    kids_real* Xnj;   ///< [N * F * F] used in Stochastic Schrodinger Equation Methods (alias Qmat)

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

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // SystemBath_H

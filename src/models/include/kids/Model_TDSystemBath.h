#ifndef Model_TDSystemBath_H
#define Model_TDSystemBath_H

#include "kids/Model.h"
#include "kids/Model_HarmonicBath.h"
#include "kids/Model_SystemBath.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(TDSystemPolicy,  //
              TD1,             // 
              TD2,             // 
              Read);           //

class Model_TDSystemBath final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_TDSystemBath() {
        appendChild(std::shared_ptr<Model_HarmonicBath>(new Model_HarmonicBath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of baths
    int Nb;     // no. of discrete modes for a bath
    int L;      // no. of nonzero variables in each interaction matrix Q

    kids_real perx, pery;
    kids_real freqd;

    // system & coupling
    span<kids_real> Hsys;  ///< [F * F] (electonic) system Hamiltonian matrix
    span<kids_real> Kmat;  ///< [N * N] (nuclear) oscillation strength matrix
    span<kids_real> Qmat;  ///< [N * F * F] coupling matrix
    span<kids_real> Q;     ///< [nbath * FF] interaction matrix for different bath
    span<kids_real> CL;    ///< [L * Nb] discretized coefficients times nonzero terms of Qj
    span<kids_real> QL;    ///< [L * Nb * F * F] save coulping matrix, each and L no. of nonzero elements
    span<kids_real> Xnj;   ///< [N * F * F] used in Stochastic Schrodinger Equation Methods (alias Qmat)

    // bath
    span<kids_real> omegas;  ///< save discrete frequencies (only for simple model, L=1)
    span<kids_real> coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    span<kids_real> x_sigma;
    span<kids_real> p_sigma;

    // integrator
    span<kids_real> x, p, m;
    span<kids_real> t;

    // model
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess;
    span<kids_real> V, dV, ddV;

    // options
    SystemPolicy::_type   system_type;
    CouplingPolicy::_type coupling_type;

    bool is_et_transform;

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_TDSystemBath_H

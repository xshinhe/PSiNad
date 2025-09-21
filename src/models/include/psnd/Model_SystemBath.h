#ifndef Model_SystemBath_H
#define Model_SystemBath_H

#include "psnd/Model.h"
#include "psnd/Model_HarmonicBath.h"
#include "psnd/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(SystemPolicy,  //
              SB,            // Spin-Boson Model:
              FMO,           // FMO model:
              SF3a,          // Singlet-Fussion Models:
              SF3b,          //
              SF3c,          //
              SF5a,          //
              SF5b,          //
              FCP,           // FCP
              AGG,           // for: B850 model
              CYC,           // for: B850 model
              Read);         //

DEFINE_POLICY(CouplingPolicy,  //
              SB,              //
              SE,              //
              Read);           //

class Model_SystemBath final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_SystemBath() {
        appendChild(std::shared_ptr<Model_HarmonicBath>(new Model_HarmonicBath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of baths
    int Nb;     // no. of discrete modes for a bath
    int L;      // no. of nonzero variables in each interaction matrix Q

    // system & coupling
    span<psnd_real> Hsys;  ///< [F * F] (electonic) system Hamiltonian matrix
    span<psnd_real> Kmat;  ///< [N * N] (nuclear) oscillation strength matrix
    span<psnd_real> Qmat;  ///< [N * F * F] coupling matrix
    span<psnd_real> Q;     ///< [nbath * FF] interaction matrix for different bath
    span<psnd_real> CL;    ///< [L * Nb] discretized coefficients times nonzero terms of Qj
    span<psnd_real> QL;    ///< [L * Nb * F * F] save coulping matrix, each and L no. of nonzero elements
    span<psnd_real> Xnj;   ///< [N * F * F] used in Stochastic Schrodinger Equation Methods (alias Qmat)

    // bath
    span<psnd_real> omegas;  ///< save discrete frequencies (only for simple model, L=1)
    span<psnd_real> coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    span<psnd_real> x_sigma;
    span<psnd_real> p_sigma;

    // integrator
    span<psnd_real> x, p, m;

    // model
    span<psnd_real> mass;
    span<psnd_real> vpes, grad, hess;
    span<psnd_real> V, dV, ddV;

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

#endif  // Model_SystemBath_H

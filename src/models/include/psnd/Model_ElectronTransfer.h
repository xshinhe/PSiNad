#ifndef Model_ElectronTransfer_H
#define Model_ElectronTransfer_H

#include "psnd/Model.h"
#include "psnd/Model_HarmonicBath.h"
#include "psnd/Policy.h"

namespace PROJECT_NS {

class Model_ElectronTransfer final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_ElectronTransfer() {
        appendChild(std::shared_ptr<Model_HarmonicBath>(new Model_HarmonicBath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of bath
    int Nb;     // discrete no.

    // system & coupling
    span<psnd_real> Hsys; /* Hamiltonian for system part */
    span<psnd_real> Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */


    psnd_real omega0;
    psnd_real lambda0;
    psnd_real coeff0;
    psnd_real beta;
    int       scan_flag;

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

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_ElectronTransfer_H

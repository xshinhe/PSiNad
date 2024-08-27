#ifndef Model_ElectronTransfer_H
#define Model_ElectronTransfer_H

#include "kids/Model.h"
#include "kids/Model_HarmonicBath.h"
#include "kids/Policy.h"

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
    span<kids_real> Hsys; /* Hamiltonian for system part */
    span<kids_real> Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */


    kids_real omega0;
    kids_real lambda0;
    kids_real coeff0;
    kids_real beta;
    int       scan_flag;

    // bath
    span<kids_real> omegas;  ///< save discrete frequencies (only for simple model, L=1)
    span<kids_real> coeffs;  ///< save coupling coefficients (only for simple model, L=1)
    span<kids_real> x_sigma;
    span<kids_real> p_sigma;

    // integrator
    span<kids_real> x, p, m;

    // model
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess;
    span<kids_real> V, dV, ddV;

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_ElectronTransfer_H

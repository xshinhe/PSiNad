#ifndef Model_ElectronTransfer_H
#define Model_ElectronTransfer_H

#include "kids/Model.h"
#include "kids/Model_Bath.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

class Model_ElectronTransfer final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_ElectronTransfer() {
        appendChild(std::shared_ptr<Model_Bath>(new Model_Bath()));  //
    }

   private:
    // parameters
    int nbath;  // no. of bath
    int Nb;     // discrete no.

    // system & coupling
    kids_real* Hsys; /* Hamiltonian for system part */
    kids_real* Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */


    kids_real omega0;
    kids_real lambda0;
    kids_real coeff0;
    kids_real beta;
    int       scan_flag;

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

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_ElectronTransfer_H

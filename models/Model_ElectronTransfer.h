#ifndef Model_ElectronTransfer_H
#define Model_ElectronTransfer_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Model_Bath.h"

namespace kids {

class Model_ElectronTransfer final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_ElectronTransfer"; }

    Model_ElectronTransfer() {
        push(std::shared_ptr<Model_Bath>(new Model_Bath()));  //
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
    int scan_flag;

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

    virtual void read_param_impl(Param* PM);
    virtual void init_data_impl(DataSet* DS);
    virtual void init_calc_impl(int stat = -1);
    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids

#endif  // Model_ElectronTransfer_H

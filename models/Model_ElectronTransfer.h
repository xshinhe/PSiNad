#ifndef Model_ElectronTransfer_H
#define Model_ElectronTransfer_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Model_Bath.h"

namespace PROJECT_NS {

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
    num_real* Hsys; /* Hamiltonian for system part */
    num_real* Q;    /* system part in interaction with different bath  [size: NvalinQ * nbath * FF] */


    num_real omega0;
    num_real lambda0;
    num_real coeff0;
    num_real beta;
    int scan_flag;

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

    virtual void read_param_impl(Param* P);
    virtual void init_data_impl(DataSet* DS);
    virtual void init_calc_impl(int stat = -1);
    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Model_ElectronTransfer_H

#ifndef MODEL_NAD1D_H
#define MODEL_NAD1D_H

#include "../core/Kernel.h"
#include "../core/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NAD1DPolicy,
              SAC,      // Tully's Single Avoid Crossing Model
              SAC2,     // Tully's Single Avoid Crossing Model (with slight revision)
              DAC,      // Tully's Doubly Avoid Crossing Model
              ECR,      // Tully's Extend Coupling Region Model
              DBG,      // (double)
              DAG,      // (double)
              DRN,      // (double) ECR
              NA_I,     // Na + I collision model
              MORSE3A,  // 3-state MORSE3A model
              MORSE3B,  // 3-state MORSE3B model
              MORSE3C,  // 3-state MORSE3C model
              MORSE15,  // 15-state MORSE model
              CL1D,     // 1-mode Caldeira-Leggett model
              JC1D,     // 1-mode Jaynes-Cummings model
              IVP1,     // iverted potential model 1
              IVP2,     // iverted potential model 2
              IVP3,     // iverted potential model 3
              IVP4);    // iverted potential model 4


class Model_NAD1D final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_NAD1D"; }

   private:
    NAD1DPolicy::_type nad1d_type;

    num_real* Hsys;

    num_real* x0;
    num_real* p0;
    num_real* x_sigma;
    num_real* p_sigma;

    // integrator
    num_real *x, *p;
    num_complex* p_sign;

    // model
    num_real* mass;
    num_real *vpes, *grad, *hess;
    num_real *V, *dV, *ddV;
    num_real* pm;

    virtual void read_param_impl(Param* P);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // MODEL_NAD1D_H

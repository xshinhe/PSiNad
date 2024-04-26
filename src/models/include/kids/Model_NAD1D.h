#ifndef MODEL_NAD1D_H
#define MODEL_NAD1D_H

#include "kids/Kernel.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NAD1DPolicy,
              PURE,      // Pure Electronic Dynamic
              SAC,       // Tully's Single Avoid Crossing Model
              SAC2,      // Tully's Single Avoid Crossing Model (with slight revision)
              SAC3,      // asymmertical Tully's Single Avoid Crossing Model
              DAC,       // Tully's Doubly Avoid Crossing Model
              ECR,       // Tully's Extend Coupling Region Model
              DBG,       // (double)
              DAG,       // (double)
              DRN,       // (double) ECR
              NA_I,      // Na + I collision model
              MORSE3A,   // 3-state MORSE3A model
              MORSE3B,   // 3-state MORSE3B model
              MORSE3C,   // 3-state MORSE3C model
              MORSE15,   // 15-state MORSE modelc
              MORSE15C,  // 15-state MORSE modelc
              MORSE15E,  // 15-state MORSE modelc
              CL1D,      // 1-mode Caldeira-Leggett model
              JC1D,      // 1-mode Jaynes-Cummings model
              IVP1,      // iverted potential model 1
              IVP2,      // iverted potential model 2
              IVP3,      // iverted potential model 3
              IVP4);     // iverted potential model 4

class Model_NAD1D final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    NAD1DPolicy::_type nad1d_type;

    kids_real* Hsys;

    kids_real* x0;
    kids_real* p0;
    kids_real* x_sigma;
    kids_real* p_sigma;

    // integrator
    kids_real *   x, *p;
    kids_complex* p_sign;

    // model
    kids_real* mass;
    kids_real *vpes, *grad, *hess;
    kids_real *V, *dV, *ddV;
    kids_real* pm;

    virtual void setInputParam_impl(std::shared_ptr<Param>& PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};


};  // namespace PROJECT_NS

#endif  // MODEL_NAD1D_H

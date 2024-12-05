#ifndef MODEL_NAD1D_H
#define MODEL_NAD1D_H

#include "kids/Model.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(NAD1DPolicy,
              PURE,      // Pure Electronic Dynamic
              SAC,       // Tully's Single Avoid Crossing Model
              SAC2,      // Tully's Single Avoid Crossing Model (with slight revision)
              SAC3,      // asymmertical Tully's Single Avoid Crossing Model
              SACX,      // another Single Avoid Crossing Model
              SACP,      // another Single Avoid Crossing Model
              DAC,       // Tully's Doubly Avoid Crossing Model
              ECR,       // Tully's Extend Coupling Region Model
              DBG,       // (double)
              DAG,       // (double)
              DRN,       // (double) ECR
              DPES,      //
              NA_I,      // Na + I collision model
              MORSE3A,   // 3-state MORSE3A model
              MORSE3B,   // 3-state MORSE3B model
              MORSE3C,   // 3-state MORSE3C model
              MORSE15,   // 15-state MORSE modelc
              MORSE15C,  // 15-state MORSE modelc
              MORSE15E,  // 15-state MORSE modelc
              CL1D,      // 1-mode Caldeira-Leggett model
              JC1D,      // 1-mode Jaynes-Cummings model
              RABI,      // 1-mode Rabi model
              IVP1,      // iverted potential model 1
              IVP2,      // iverted potential model 2
              IVP3,      // iverted potential model 3
              IVP4);     // iverted potential model 4

class Model_NAD1D final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    NAD1DPolicy::_type nad1d_type;

    span<kids_real> Hsys;

    span<kids_real> x0;
    span<kids_real> p0;
    span<kids_real> x_sigma;
    span<kids_real> p_sigma;

    // integrator
    span<kids_real>    x, p;
    span<kids_complex> p_sign;

    // model
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess;
    span<kids_real> V, dV, ddV;
    span<kids_real> pm;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};


};  // namespace PROJECT_NS

#endif  // MODEL_NAD1D_H

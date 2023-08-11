#ifndef Kernel_NADForce_H
#define Kernel_NADForce_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Representation.h"

namespace PROJECT_NS {

DEFINE_POLICY(NADForcePolicy,  //
              BO,              //
              EHR,             //
              ELSE);           //

class Kernel_NADForce : public Kernel {
   public:
    static NADForcePolicy::_type NADForce_type;

    inline virtual const std::string name() { return "Kernel_NADForce"; }

   private:
    num_real *f, *grad, *dV, *dE, *Force, *T;

    bool BATH_FORCE_OPT;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_NADForce_H

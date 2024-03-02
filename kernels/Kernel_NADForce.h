#ifndef Kernel_NADForce_H
#define Kernel_NADForce_H

#include "../core/Kernel.h"
#include "../core/Policy.h"

namespace kids {

DEFINE_POLICY(NADForcePolicy,  //
              EHR,             //
              BO,              //
              CV,              //
              BOSD,            //
              CVSD,            //
              ELSE);           //

namespace FORCE_OPT {
extern bool BATH_FORCE_BILINEAR;
extern int nbath;
extern int Nb;
};  // namespace FORCE_OPT

class Kernel_NADForce : public Kernel {
   public:
    static NADForcePolicy::_type NADForce_type;
    inline virtual const std::string name() { return "Kernel_NADForce"; }

   private:
    bool offd_projected;

    kids_real *f, *grad, *dV, *dE, *Force, *T;
    kids_real *p, *m;
    kids_real *fadd, *fproj;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace kids

#endif  // Kernel_NADForce_H

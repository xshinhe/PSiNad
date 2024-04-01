#ifndef Model_LVCM_H
#define Model_LVCM_H

#include "../core/Kernel.h"
#include "../core/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(LVCMPolicy,  //
              PYR3,        //
              PYR3_SPEC,   //
              PYR4,        //
              PYR4_SPEC,   //
              PYR24,       //
              CRC2,        //
              CRC5,        //
              BEN5,        //
              BUTA5,       //
              PENTA5,      //
              CED2,        //
              CED3,        //
              PYR2CED,     //
              Read);       //

class Model_LVCM final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_LVCM"; }

    Model_LVCM(){};

   private:
    LVCMPolicy::_type lvcm_type;
    int               N_coup;
    int               N_mode;
    bool              classical_bath;


    kids_real* Hsys;
    kids_real *kcoeff, *lcoeff;

    kids_real* x_sigma;
    kids_real* p_sigma;
    kids_real* x0;
    kids_real* p0;

    // integrator
    kids_real *x, *p, *m, *w;

    // model
    kids_real* mass;
    kids_real *vpes, *grad, *hess;
    kids_real *V, *dV, *ddV;

    // int N_ligh;
    // N = N_mode + N_coup + N_ligh

    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat);

    int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Model_LVCM_H

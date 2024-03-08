#ifndef Kernel_Update_T_H
#define Kernel_Update_T_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_T : public Kernel {
   public:
    Kernel_Update_T(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_T"; }

   private:
    double *p, *m;
    // for Langevin
    double *c1, *c2p;
    // for NHC
    double *nhc_x, *nhc_p, *nhc_G, *nhc_Q;
    //
    double scale, *dt_ptr;
    double beta;
    double gammal;
    double randu;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Update_T_H

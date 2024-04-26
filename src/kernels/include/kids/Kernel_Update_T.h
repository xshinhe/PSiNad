#ifndef Kernel_Update_T_H
#define Kernel_Update_T_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_T : public Kernel {
   public:
    Kernel_Update_T(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

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

    virtual void setInputParam_impl(std::shared_ptr<Param> &PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> &DS);

    virtual Status &executeKernel_impl(Status &stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Update_T_H

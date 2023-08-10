#ifndef Kernel_Representation_H
#define Kernel_Representation_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "../core/linalg.h"

namespace PROJECT_NS {

DEFINE_POLICY(RepresentationPolicy,
              Diabatic,   // diabatic representation
              Adiabatic,  // adiabtic representation
              Onthefly,   // onthefly adiabtic representation
              Force,      // adiabtic representation
              Density     // adiabtic representation
);

/**
 * @brief      Kernel_Representation solves basis transformation for electronic problems
 * the data is storaged on field "sol.rep." (might be an alias)
 */
class Kernel_Representation final : public Kernel {
   public:
    static RepresentationPolicy::_type representation_type;

    inline virtual const std::string name() { return "Kernel_Representation"; }

   private:
    // int Nc, N, F, FF, NFF, NNFF;
    std::complex<double>* c;
    std::complex<double>* rho;

    double *V, *dV, *ddV;
    double *E, *T, *dE, *ddE;
    double* L;
    std::complex<double>*R, *dL, *ddL;
    std::complex<double>*H, *dH, *ddH;

    double *x, *p, *m;
    double *Told, *Tnew, *TtTold, *Enew;
    double* Matr_tmp;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};



};  // namespace PROJECT_NS

#endif  // Kernel_Representation_H

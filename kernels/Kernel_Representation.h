#ifndef Kernel_Representation_H
#define Kernel_Representation_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "../core/linalg.h"

namespace PROJECT_NS {

DEFINE_POLICY(RepresentationPolicy,
              Diabatic,   // diabatic representation
              Adiabatic,  // adiabtic representation
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
    static RepresentationPolicy::_type ini_repr_type;
    static RepresentationPolicy::_type ele_repr_type;
    static RepresentationPolicy::_type nuc_repr_type;
    static RepresentationPolicy::_type tcf_repr_type;
    static bool onthefly;

    inline virtual const std::string name() { return "Kernel_Representation"; }

   private:
    bool do_refer;
    bool phase_correction;

    // int Nc, N, F, FF, NFF, NNFF;
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

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};



};  // namespace PROJECT_NS

#endif  // Kernel_Representation_H

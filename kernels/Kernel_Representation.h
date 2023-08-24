#ifndef Kernel_Representation_H
#define Kernel_Representation_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "../core/linalg.h"

namespace PROJECT_NS {

DEFINE_POLICY(RepresentationPolicy,
              Diabatic,   // diabatic representation
              Adiabatic,  // adiabtic representation
              Force,      // @developing in future
              Density     // @developing in future
);

/**
 * @brief  Kernel_Representationfor solving basis transformation for electronic problems
 *
 */
class Kernel_Representation final : public Kernel {
   public:
    static RepresentationPolicy::_type representation_type;  ///< root representation
    static RepresentationPolicy::_type inp_repr_type;        ///< representatioin for input (quantities)
    static RepresentationPolicy::_type ele_repr_type;        ///< representation for electronic dynamics
    static RepresentationPolicy::_type nuc_repr_type;        ///< representation for nuclear dynamics
    static RepresentationPolicy::_type tcf_repr_type;        ///< representation for intrinsic calculation of TCF
    static bool onthefly;                                    ///< flag indicated the ab inition calculation

    inline virtual const std::string name() { return "Kernel_Representation"; }

   private:
    bool do_refer;
    bool phase_correction;

    double *V, *dV, *ddV;
    double *E, *T, *Told, *dE, *ddE;
    double *L;
    std::complex<double> *R, *dL, *ddL;
    std::complex<double> *H, *dH, *ddH;

    double *x, *p, *m;
    int *occ_nuc;
    num_complex *rho_ele;
    double *ve, *vedE, *TtTold;

    virtual void read_param_impl(Param *PM);

    virtual void init_data_impl(DataSet *DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Representation_H

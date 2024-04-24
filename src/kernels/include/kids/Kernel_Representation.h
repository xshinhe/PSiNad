#ifndef Kernel_Representation_H
#define Kernel_Representation_H

#include "kids/Kernel.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(RepresentationPolicy,
              Diabatic,   // diabatic representation
              Adiabatic,  // adiabtic representation
              Force,      // @developing in future
              Density     // @developing in future
);

DEFINE_POLICY(SpacePolicy,
              H,  // Hilbert space
              L   // Liouvillian space
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
    static bool                        onthefly;             ///< flag indicated the ab inition calculation

    inline virtual const std::string name() { return "Kernel_Representation"; }

    static int transform(kids_complex *A, kids_real *T, int fdim,  //
                         RepresentationPolicy::_type from, RepresentationPolicy::_type to, SpacePolicy::_type Stype);

   private:
    bool do_refer;
    bool phase_correction;
    bool basis_switch;

    double *              V, *dV, *ddV;
    double *              E, *T, *Told, *dE, *ddE;
    double *              L;
    std::complex<double> *R, *dL, *ddL;
    std::complex<double> *H, *dH, *ddH;

    double *E_copy;

    double *      x, *p, *m;
    int *         occ_nuc;
    kids_complex *rho_ele;
    double *      ve, *vedE, *TtTold;

    virtual void setInputParam_impl(Param *PM);

    virtual void setInputDataSet_impl(DataSet *DS);

    virtual Status &initializeKernel_impl(Status &stat);

    virtual Status &executeKernel_impl(Status &stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Representation_H

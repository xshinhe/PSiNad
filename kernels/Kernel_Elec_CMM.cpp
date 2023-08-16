#include "Kernel_Elec_CMM.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Elec.h"
#include "Kernel_Random.h"


#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

namespace PROJECT_NS {

/**
 * @brief function for gamma_wigner: gamma = (sqrt(1+F)-1)/F
 */
double Kernel_Elec_CMM::gamma_wigner(int Fdim) { return (sqrt((num_real) (Fdim) + 1) - 1) / Fdim; }

/**
 * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
 */
int Kernel_Elec_CMM::c_sphere(num_complex* c, int fdim) {
    num_real rpart, ipart;
    for (int i = 0; i < fdim; ++i) {
        Kernel_Random::rand_gaussian(&rpart);
        Kernel_Random::rand_gaussian(&ipart);
        c[i] = num_complex(rpart, ipart);
    }
    num_complex norm2;
    ARRAY_MATMUL_TRANS1(&norm2, c, c, 1, fdim, 1);
    double l = std::sqrt(std::real(norm2));
    for (int i = 0; i < fdim; ++i) c[i] /= l;
    return 0;
}

void Kernel_Elec_CMM::read_param_impl(Param* PM) {
    gamma1  = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));
    gamma2  = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1     = (1 + Dimension::F * gamma1);
    xi2     = (1 + Dimension::F * gamma2);
    use_cv  = PM->get<bool>("use_cv", LOC(), false);
    use_wmm = PM->get<bool>("use_wmm", LOC(), false);  // @disable
}

void Kernel_Elec_CMM::init_calc_impl(int stat) {
    Kernel_Elec::w[0] = num_complex(Dimension::F);
    c_sphere(Kernel_Elec::c, Dimension::F);

    *Kernel_Elec::occ_nuc = Kernel_Elec::occ0;                                          // useless
    Kernel_Elec::ker_from_c(Kernel_Elec::rho_ele, Kernel_Elec::c, 1, 0, Dimension::F);  // single-rank

    Kernel_Elec::ker_from_rho(Kernel_Elec::rho_nuc, Kernel_Elec::rho_ele, xi1, gamma1, Dimension::F);
    if (use_cv) {
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            Kernel_Elec::rho_nuc[ii] = (i == Kernel_Elec::occ0) ? phys::math::iu : phys::math::iz;
        }
    }

    exec_kernel(stat);
}

int Kernel_Elec_CMM::exec_kernel_impl(int stat) {
    // calc TCF kernels
    Kernel_Elec::ker_from_rho(Kernel_Elec::K1, Kernel_Elec::rho_ele, xi1, gamma1, Dimension::F);
    Kernel_Elec::ker_from_rho(Kernel_Elec::K2, Kernel_Elec::rho_ele, xi2, gamma2, Dimension::F);
    return 0;
}

};  // namespace PROJECT_NS

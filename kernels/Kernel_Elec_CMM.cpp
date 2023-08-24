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


double Kernel_Elec_CMM::gamma_wigner(int Fdim) { return (sqrt((num_real) (Fdim) + 1) - 1) / Fdim; }

/**
 * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
 */
int Kernel_Elec_CMM::c_sphere(num_complex* c, int fdim) {
    Kernel_Random::rand_gaussian(reinterpret_cast<double*>(c), 2 * fdim);
    double l = std::sqrt(std::real(ARRAY_INNER_CONJ1(c, c, fdim)));
    for (int i = 0; i < fdim; ++i) c[i] /= l;
    return 0;
}

void Kernel_Elec_CMM::read_param_impl(Param* PM) {
    gamma1  = PM->get<num_real>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(Dimension::F));
    gamma2  = (1 - gamma1) / (1.0f + Dimension::F * gamma1);
    xi1     = (1 + Dimension::F * gamma1);
    xi2     = (1 + Dimension::F * gamma2);
    use_cv  = PM->get<bool>("use_cv", LOC(), false);
    use_wmm = PM->get<bool>("use_wmm", LOC(), false);
}

void Kernel_Elec_CMM::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* w       = Kernel_Elec::w + iP;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc         = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0]     = num_complex(Dimension::F);                     ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;                             ///< initial occupation
        c_sphere(c, Dimension::F);                                ///< initial c on standard sphere
        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele
        Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, xi1, gamma1, Dimension::F, use_cv, *occ_nuc);  ///< initial rho_nuc
        ARRAY_EYE(U, Dimension::F);  ///< initial propagator
    }

    Kernel_Elec::c_init       = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->set("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::rho_nuc_init = _DataSet->set("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    exec_kernel(stat);
}

int Kernel_Elec_CMM::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* c            = Kernel_Elec::c + iP * Dimension::F;
        num_complex* c_init       = Kernel_Elec::c_init + iP * Dimension::F;
        num_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        num_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* rho_nuc_init = Kernel_Elec::rho_nuc_init + iP * Dimension::FF;
        num_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        num_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;

        ARRAY_MATMUL(c, U, c_init, Dimension::F, Dimension::F, 1);
        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele_init, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc_init, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);

        Kernel_Elec::ker_from_rho(K1, rho_ele, xi1, gamma1, Dimension::F);
        Kernel_Elec::ker_from_rho(K2, rho_ele, xi2, gamma2, Dimension::F);
    }
    return 0;
}

};  // namespace PROJECT_NS

#include "Kernel_Elec.h"

#include "../core/linalg.h"
#include "../core/vars_list.h"
#include "Kernel_Declare.h"
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

void Kernel_Elec::read_param_impl(Param* PM) {
    occ0 = PM->get<int>("occ", LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Dimension::F) throw std::runtime_error("occ >= F");
}

void Kernel_Elec::init_data_impl(DataSet* DS) {
    U        = DS->def(DATA::integrator::U);
    c        = DS->def(DATA::integrator::c);
    rho_ele  = DS->def(DATA::integrator::rho_ele);
    rho_dual = DS->def(DATA::integrator::rho_dual);
    T        = DS->def(DATA::model::rep::T);
    occ_nuc  = DS->def(DATA::integrator::occ_nuc);
    rho_nuc  = DS->def(DATA::integrator::rho_nuc);

    w    = DS->def(DATA::integrator::w);
    w_CC = DS->def(DATA::integrator::w_CC);
    w_CP = DS->def(DATA::integrator::w_CP);
    w_PP = DS->def(DATA::integrator::w_PP);
    w_AA = DS->def(DATA::integrator::w_AA);
    w_AD = DS->def(DATA::integrator::w_AD);
    w_DD = DS->def(DATA::integrator::w_DD);
    wz_A = DS->def(DATA::integrator::wz_A);
    wz_D = DS->def(DATA::integrator::wz_D);
    ww_A = DS->def(DATA::integrator::ww_A);
    ww_D = DS->def(DATA::integrator::ww_D);

    K0   = DS->def(DATA::integrator::K0);
    K1   = DS->def(DATA::integrator::K1);
    K2   = DS->def(DATA::integrator::K2);
    K1QA = DS->def(DATA::integrator::K1QA);
    K2QA = DS->def(DATA::integrator::K2QA);
    K1DA = DS->def(DATA::integrator::K1DA);
    K2DA = DS->def(DATA::integrator::K2DA);
    K1QD = DS->def(DATA::integrator::K1QD);
    K2QD = DS->def(DATA::integrator::K2QD);
    K1DD = DS->def(DATA::integrator::K1DD);
    K2DD = DS->def(DATA::integrator::K2DD);

    // read input operator
    OpA = DS->def(DATA::integrator::OpA);
    OpB = DS->def(DATA::integrator::OpB);
}

void Kernel_Elec::init_calc_impl(int stat) {
    exec_kernel_impl(stat);

    double unit = 1.0e0;
    _DataSet->def("integrator.1", &unit);
    _DataSet->def("init.1", &unit);
    _DataSet->def("init.w", w, Dimension::P);
    _DataSet->def("init.K0", K0, Dimension::PFF);
    _DataSet->def("init.K1", K1, Dimension::PFF);
    _DataSet->def("init.K2", K2, Dimension::PFF);
    _DataSet->def("init.K1QA", K1QA, Dimension::PFF);
    _DataSet->def("init.K2QA", K2QA, Dimension::PFF);
    _DataSet->def("init.K1DA", K1DA, Dimension::PFF);
    _DataSet->def("init.K2DA", K2DA, Dimension::PFF);
    _DataSet->def("init.K1QD", K1QD, Dimension::PFF);
    _DataSet->def("init.K2QD", K2QD, Dimension::PFF);
    _DataSet->def("init.K1DD", K1DD, Dimension::PFF);
    _DataSet->def("init.K2DD", K2DD, Dimension::PFF);
}

int Kernel_Elec::exec_kernel_impl(int stat) { return 0; }

/**
 * @brief convert c (electonic amplititude) to kernel (affine map of the density)
 */
int Kernel_Elec::ker_from_c(kids_complex* ker, kids_complex* c, kids_real xi, kids_real gamma, int fdim) {
    ARRAY_OUTER_TRANS2(ker, c, c, fdim, fdim);
    for (int i = 0, idx = 0; i < fdim; ++i) {
        for (int j = 0; j < fdim; ++j, ++idx) {
            ker[idx] *= xi;
            if (i == j) ker[idx] -= phys::math::iu * gamma;
        }
    }
    return 0;
}


/**
 * @brief convert c (electonic amplititude) to kernel (affine map of the density)
 */
int Kernel_Elec::ker_from_rho(kids_complex* ker, kids_complex* rho, kids_real xi, kids_real gamma, int fdim,
                              bool quantize, int occ) {
    for (int i = 0, ij = 0; i < fdim; ++i) {
        for (int j = 0; j < fdim; ++j, ++ij) {
            ker[ij] = xi * rho[ij];
            if (i == j) ker[ij] = (quantize) ? (i == occ ? 1.0e0 : 0.0e0) : (ker[ij] - gamma);
        }
    }
    return 0;
}

int Kernel_Elec::occ0;

// two densities for dynamic
kids_complex* Kernel_Elec::U;
kids_complex* Kernel_Elec::c;
kids_complex* Kernel_Elec::c_init;
kids_complex* Kernel_Elec::rho_ele;
kids_complex* Kernel_Elec::rho_ele_init;
kids_complex* Kernel_Elec::rho_dual;
kids_complex* Kernel_Elec::rho_dual_init;
kids_real*    Kernel_Elec::T;
kids_real*    Kernel_Elec::T_init;

int*          Kernel_Elec::occ_nuc;
kids_complex* Kernel_Elec::rho_nuc;
kids_complex* Kernel_Elec::rho_nuc_init;

// time correlation function
kids_complex* Kernel_Elec::w;
kids_complex* Kernel_Elec::w_CC;
kids_complex* Kernel_Elec::w_CP;
kids_complex* Kernel_Elec::w_PP;
kids_complex* Kernel_Elec::w_AA;
kids_complex* Kernel_Elec::w_AD;
kids_complex* Kernel_Elec::w_DD;
kids_complex* Kernel_Elec::wz_A;
kids_complex* Kernel_Elec::wz_D;
kids_complex* Kernel_Elec::ww_A;
kids_complex* Kernel_Elec::ww_D;
kids_complex* Kernel_Elec::ww_A_init;
kids_complex* Kernel_Elec::ww_D_init;
kids_complex* Kernel_Elec::K0;
kids_complex* Kernel_Elec::K1;
kids_complex* Kernel_Elec::K2;
kids_complex* Kernel_Elec::K1QA;
kids_complex* Kernel_Elec::K2QA;
kids_complex* Kernel_Elec::K1DA;
kids_complex* Kernel_Elec::K2DA;
kids_complex* Kernel_Elec::K1QD;
kids_complex* Kernel_Elec::K2QD;
kids_complex* Kernel_Elec::K1DD;
kids_complex* Kernel_Elec::K2DD;

kids_complex* Kernel_Elec::OpA;
kids_complex* Kernel_Elec::OpB;
kids_complex* Kernel_Elec::TrK1A;
kids_complex* Kernel_Elec::TrK2B;

};  // namespace PROJECT_NS

#include "Kernel_Elec.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Random.h"

namespace PROJECT_NS {

void Kernel_Elec::read_param_impl(Param* PM) {
    occ0 = PM->get<int>("occ", LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Dimension::F) throw std::runtime_error("occ >= F");
}

void Kernel_Elec::init_data_impl(DataSet* DS) {
    U        = DS->reg<num_complex>("integrator.U", Dimension::PFF);
    c        = DS->reg<num_complex>("integrator.c", Dimension::PF);
    rho_ele  = DS->reg<num_complex>("integrator.rho_ele", Dimension::PFF);
    rho_dual = DS->reg<num_complex>("integrator.rho_dual", Dimension::PFF);
    T        = DS->reg<num_real>("model.rep.T", Dimension::PFF);

    occ_nuc = DS->reg<int>("integrator.occ_nuc", Dimension::P);
    rho_nuc = DS->reg<num_complex>("integrator.rho_nuc", Dimension::PFF);

    w    = DS->reg<num_complex>("integrator.w", Dimension::P);
    w_CC = DS->reg<num_complex>("integrator.w_CC", Dimension::P);
    w_CP = DS->reg<num_complex>("integrator.w_CP", Dimension::P);
    w_PP = DS->reg<num_complex>("integrator.w_PP", Dimension::P);
    w_AA = DS->reg<num_complex>("integrator.w_AA", Dimension::P);
    w_AD = DS->reg<num_complex>("integrator.w_AD", Dimension::P);
    w_DD = DS->reg<num_complex>("integrator.w_DD", Dimension::P);
    wz_A = DS->reg<num_complex>("integrator.wz_A", Dimension::P);
    wz_D = DS->reg<num_complex>("integrator.wz_D", Dimension::P);
    ww_A = DS->reg<num_complex>("integrator.ww_A", Dimension::P);
    ww_D = DS->reg<num_complex>("integrator.ww_D", Dimension::P);

    K0    = DS->reg<num_complex>("integrator.K0", Dimension::PFF);
    K0dia = DS->reg<num_complex>("integrator.K0dia", Dimension::PF);
    K0occ = DS->reg<num_complex>("integrator.K0occ", Dimension::P);

    K1    = DS->reg<num_complex>("integrator.K1", Dimension::PFF);
    K1dia = DS->reg<num_complex>("integrator.K1dia", Dimension::PF);
    K1occ = DS->reg<num_complex>("integrator.K1occ", Dimension::P);

    K2    = DS->reg<num_complex>("integrator.K2", Dimension::PFF);
    K2dia = DS->reg<num_complex>("integrator.K2dia", Dimension::PF);
    K2occ = DS->reg<num_complex>("integrator.K2occ", Dimension::P);

    K1QA    = DS->reg<num_complex>("integrator.K1QA", Dimension::PFF);
    K1QAdia = DS->reg<num_complex>("integrator.K1QAdia", Dimension::PF);
    K1QAocc = DS->reg<num_complex>("integrator.K1QAocc", Dimension::P);
    K2QA    = DS->reg<num_complex>("integrator.K2QA", Dimension::PFF);
    K2QAdia = DS->reg<num_complex>("integrator.K2QAdia", Dimension::PF);
    K2QAocc = DS->reg<num_complex>("integrator.K2QAocc", Dimension::P);
    K1DA    = DS->reg<num_complex>("integrator.K1DA", Dimension::PFF);
    K1DAdia = DS->reg<num_complex>("integrator.K1DAdia", Dimension::PF);
    K1DAocc = DS->reg<num_complex>("integrator.K1DAocc", Dimension::P);
    K2DA    = DS->reg<num_complex>("integrator.K2DA", Dimension::PFF);
    K2DAdia = DS->reg<num_complex>("integrator.K2DAdia", Dimension::PF);
    K2DAocc = DS->reg<num_complex>("integrator.K2DAocc", Dimension::P);

    K1QD    = DS->reg<num_complex>("integrator.K1QD", Dimension::PFF);
    K1QDdia = DS->reg<num_complex>("integrator.K1QDdia", Dimension::PF);
    K1QDocc = DS->reg<num_complex>("integrator.K1QDocc", Dimension::P);
    K2QD    = DS->reg<num_complex>("integrator.K2QD", Dimension::PFF);
    K2QDdia = DS->reg<num_complex>("integrator.K2QDdia", Dimension::PF);
    K2QDocc = DS->reg<num_complex>("integrator.K2QDocc", Dimension::P);
    K1DD    = DS->reg<num_complex>("integrator.K1DD", Dimension::PFF);
    K1DDdia = DS->reg<num_complex>("integrator.K1DDdia", Dimension::PF);
    K1DDocc = DS->reg<num_complex>("integrator.K1DDocc", Dimension::P);
    K2DD    = DS->reg<num_complex>("integrator.K2DD", Dimension::PFF);
    K2DDdia = DS->reg<num_complex>("integrator.K2DDdia", Dimension::PF);
    K2DDocc = DS->reg<num_complex>("integrator.K2DDocc", Dimension::P);

    // read input operator
    OpA = DS->reg<num_complex>("integrator.OpA", Dimension::FF);
    OpB = DS->reg<num_complex>("integrator.OpB", Dimension::FF);
}

void Kernel_Elec::init_calc_impl(int stat) {
    exec_kernel_impl(stat);

    double unit = 1.0e0;
    _DataSet->set("integrator.1", &unit);
    _DataSet->set("init.1", &unit);
    _DataSet->set("init.w", w, Dimension::P);
    _DataSet->set("init.K0", K0, Dimension::PFF);
    _DataSet->set("init.K0dia", K0dia, Dimension::PF);
    _DataSet->set("init.K0occ", K0occ, Dimension::P);
    _DataSet->set("init.K1", K1, Dimension::PFF);
    _DataSet->set("init.K1dia", K1dia, Dimension::PF);
    _DataSet->set("init.K1occ", K1occ, Dimension::P);
    _DataSet->set("init.K2", K2, Dimension::PFF);
    _DataSet->set("init.K2dia", K2dia, Dimension::PF);
    _DataSet->set("init.K2occ", K2occ, Dimension::P);
    _DataSet->set("init.K1QA", K1QA, Dimension::PFF);
    _DataSet->set("init.K1QAdia", K1QAdia, Dimension::PF);
    _DataSet->set("init.K1QAocc", K1QAocc, Dimension::P);
    _DataSet->set("init.K2QA", K2QA, Dimension::PFF);
    _DataSet->set("init.K2QAdia", K2QAdia, Dimension::PF);
    _DataSet->set("init.K2QAocc", K2QAocc, Dimension::P);
    _DataSet->set("init.K1DA", K1DA, Dimension::PFF);
    _DataSet->set("init.K1DAdia", K1DAdia, Dimension::PF);
    _DataSet->set("init.K1DAocc", K1DAocc, Dimension::P);
    _DataSet->set("init.K2DA", K2DA, Dimension::PFF);
    _DataSet->set("init.K2DAdia", K2DAdia, Dimension::PF);
    _DataSet->set("init.K2DAocc", K2DAocc, Dimension::P);
    _DataSet->set("init.K1QD", K1QD, Dimension::PFF);
    _DataSet->set("init.K1QDdia", K1QDdia, Dimension::PF);
    _DataSet->set("init.K1QDocc", K1QDocc, Dimension::P);
    _DataSet->set("init.K2QD", K2QD, Dimension::PFF);
    _DataSet->set("init.K2QDdia", K2QDdia, Dimension::PF);
    _DataSet->set("init.K2QDocc", K2QDocc, Dimension::P);
    _DataSet->set("init.K1DD", K1DD, Dimension::PFF);
    _DataSet->set("init.K1DDdia", K1DDdia, Dimension::PF);
    _DataSet->set("init.K1DDocc", K1DDocc, Dimension::P);
    _DataSet->set("init.K2DD", K2DD, Dimension::PFF);
    _DataSet->set("init.K2DDdia", K2DDdia, Dimension::PF);
    _DataSet->set("init.K2DDocc", K2DDocc, Dimension::P);
}

int Kernel_Elec::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* K0occ   = this->K0occ + iP;
        num_complex* K0dia   = this->K0dia + iP * Dimension::F;
        num_complex* K0      = this->K0 + iP * Dimension::FF;
        num_complex* K1occ   = this->K1occ + iP;
        num_complex* K1dia   = this->K1dia + iP * Dimension::F;
        num_complex* K1      = this->K1 + iP * Dimension::FF;
        num_complex* K2occ   = this->K2occ + iP;
        num_complex* K2dia   = this->K2dia + iP * Dimension::F;
        num_complex* K2      = this->K2 + iP * Dimension::FF;
        num_complex* K1QAocc = this->K1QAocc + iP;
        num_complex* K1QAdia = this->K1QAdia + iP * Dimension::F;
        num_complex* K1QA    = this->K1QA + iP * Dimension::FF;
        num_complex* K2QAocc = this->K2QAocc + iP;
        num_complex* K2QAdia = this->K2QAdia + iP * Dimension::F;
        num_complex* K2QA    = this->K2QA + iP * Dimension::FF;
        num_complex* K1DAocc = this->K1DAocc + iP;
        num_complex* K1DAdia = this->K1DAdia + iP * Dimension::F;
        num_complex* K1DA    = this->K1DA + iP * Dimension::FF;
        num_complex* K2DAocc = this->K2DAocc + iP;
        num_complex* K2DAdia = this->K2DAdia + iP * Dimension::F;
        num_complex* K2DA    = this->K2DA + iP * Dimension::FF;
        num_complex* K1QDocc = this->K1QDocc + iP;
        num_complex* K1QDdia = this->K1QDdia + iP * Dimension::F;
        num_complex* K1QD    = this->K1QD + iP * Dimension::FF;
        num_complex* K2QDocc = this->K2QDocc + iP;
        num_complex* K2QDdia = this->K2QDdia + iP * Dimension::F;
        num_complex* K2QD    = this->K2QD + iP * Dimension::FF;
        num_complex* K1DDocc = this->K1DDocc + iP;
        num_complex* K1DDdia = this->K1DDdia + iP * Dimension::F;
        num_complex* K1DD    = this->K1DD + iP * Dimension::FF;
        num_complex* K2DDocc = this->K2DDocc + iP;
        num_complex* K2DDdia = this->K2DDdia + iP * Dimension::F;
        num_complex* K2DD    = this->K2DD + iP * Dimension::FF;

        ///////////////////////////////////////////////////////

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K0dia[i] = K0[ii];
        K0occ[0] = K0[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1dia[i] = K1[ii];
        K1occ[0] = K1[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2dia[i] = K2[ii];
        K2occ[0] = K2[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1QAdia[i] = K1QA[ii];
        K1QAocc[0] = K1QA[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2QAdia[i] = K2QA[ii];
        K2QAocc[0] = K2QA[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1DAdia[i] = K1DA[ii];
        K1DAocc[0] = K1DA[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2DAdia[i] = K2DA[ii];
        K2DAocc[0] = K2DA[occ0 * Dimension::Fadd1];
        //
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1QDdia[i] = K1QD[ii];
        K1QDocc[0] = K1QD[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2QDdia[i] = K2QD[ii];
        K2QDocc[0] = K2QD[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1DDdia[i] = K1DD[ii];
        K1DDocc[0] = K1DD[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2DDdia[i] = K2DD[ii];
        K2DDocc[0] = K2DD[occ0 * Dimension::Fadd1];
    }
    return 0;
}

/**
 * @brief convert c (electonic amplititude) to kernel (affine map of the density)
 */
int Kernel_Elec::ker_from_c(num_complex* ker, num_complex* c, num_real xi, num_real gamma, int fdim) {
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
int Kernel_Elec::ker_from_rho(num_complex* ker, num_complex* rho, num_real xi, num_real gamma, int fdim, bool quantize,
                              int occ) {
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
num_complex* Kernel_Elec::U;
num_complex* Kernel_Elec::c;
num_complex* Kernel_Elec::c_init;
num_complex* Kernel_Elec::rho_ele;
num_complex* Kernel_Elec::rho_ele_init;
num_complex* Kernel_Elec::rho_dual;
num_complex* Kernel_Elec::rho_dual_init;
num_real* Kernel_Elec::T;
num_real* Kernel_Elec::T_init;

int* Kernel_Elec::occ_nuc;
num_complex* Kernel_Elec::rho_nuc;
num_complex* Kernel_Elec::rho_nuc_init;

// time correlation function
num_complex* Kernel_Elec::w;
num_complex* Kernel_Elec::w_CC;
num_complex* Kernel_Elec::w_CP;
num_complex* Kernel_Elec::w_PP;
num_complex* Kernel_Elec::w_AA;
num_complex* Kernel_Elec::w_AD;
num_complex* Kernel_Elec::w_DD;
num_complex* Kernel_Elec::wz_A;
num_complex* Kernel_Elec::wz_D;
num_complex* Kernel_Elec::ww_A;
num_complex* Kernel_Elec::ww_D;
num_complex* Kernel_Elec::ww_A_init;
num_complex* Kernel_Elec::ww_D_init;
num_complex* Kernel_Elec::K0;
num_complex* Kernel_Elec::K0occ;
num_complex* Kernel_Elec::K0dia;
num_complex* Kernel_Elec::K1;
num_complex* Kernel_Elec::K1occ;
num_complex* Kernel_Elec::K1dia;
num_complex* Kernel_Elec::K2;
num_complex* Kernel_Elec::K2occ;
num_complex* Kernel_Elec::K2dia;
num_complex* Kernel_Elec::K1QA;
num_complex* Kernel_Elec::K1QAocc;
num_complex* Kernel_Elec::K1QAdia;
num_complex* Kernel_Elec::K2QA;
num_complex* Kernel_Elec::K2QAocc;
num_complex* Kernel_Elec::K2QAdia;
num_complex* Kernel_Elec::K1DA;
num_complex* Kernel_Elec::K1DAocc;
num_complex* Kernel_Elec::K1DAdia;
num_complex* Kernel_Elec::K2DA;
num_complex* Kernel_Elec::K2DAocc;
num_complex* Kernel_Elec::K2DAdia;
num_complex* Kernel_Elec::K1QD;
num_complex* Kernel_Elec::K1QDocc;
num_complex* Kernel_Elec::K1QDdia;
num_complex* Kernel_Elec::K2QD;
num_complex* Kernel_Elec::K2QDocc;
num_complex* Kernel_Elec::K2QDdia;
num_complex* Kernel_Elec::K1DD;
num_complex* Kernel_Elec::K1DDocc;
num_complex* Kernel_Elec::K1DDdia;
num_complex* Kernel_Elec::K2DD;
num_complex* Kernel_Elec::K2DDocc;
num_complex* Kernel_Elec::K2DDdia;

num_complex* Kernel_Elec::OpA;
num_complex* Kernel_Elec::OpB;
num_complex* Kernel_Elec::TrK1A;
num_complex* Kernel_Elec::TrK2B;

};  // namespace PROJECT_NS

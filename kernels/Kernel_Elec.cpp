#include "Kernel_Elec.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_Random.h"

namespace kids {

void Kernel_Elec::read_param_impl(Param* PM) {
    occ0 = PM->get<int>("occ", LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Dimension::F) throw std::runtime_error("occ >= F");
}

void Kernel_Elec::init_data_impl(DataSet* DS) {
    U        = DS->def<kids_complex>("integrator.U", Dimension::PFF);
    c        = DS->def<kids_complex>("integrator.c", Dimension::PF);
    rho_ele  = DS->def<kids_complex>("integrator.rho_ele", Dimension::PFF);
    rho_dual = DS->def<kids_complex>("integrator.rho_dual", Dimension::PFF);
    T        = DS->def<kids_real>("model.rep.T", Dimension::PFF);

    occ_nuc = DS->def<int>("integrator.occ_nuc", Dimension::P);
    rho_nuc = DS->def<kids_complex>("integrator.rho_nuc", Dimension::PFF);

    w    = DS->def<kids_complex>("integrator.w", Dimension::P);
    w_CC = DS->def<kids_complex>("integrator.w_CC", Dimension::P);
    w_CP = DS->def<kids_complex>("integrator.w_CP", Dimension::P);
    w_PP = DS->def<kids_complex>("integrator.w_PP", Dimension::P);
    w_AA = DS->def<kids_complex>("integrator.w_AA", Dimension::P);
    w_AD = DS->def<kids_complex>("integrator.w_AD", Dimension::P);
    w_DD = DS->def<kids_complex>("integrator.w_DD", Dimension::P);
    wz_A = DS->def<kids_complex>("integrator.wz_A", Dimension::P);
    wz_D = DS->def<kids_complex>("integrator.wz_D", Dimension::P);
    ww_A = DS->def<kids_complex>("integrator.ww_A", Dimension::P);
    ww_D = DS->def<kids_complex>("integrator.ww_D", Dimension::P);

    K0    = DS->def<kids_complex>("integrator.K0", Dimension::PFF);
    K0dia = DS->def<kids_complex>("integrator.K0dia", Dimension::PF);
    K0occ = DS->def<kids_complex>("integrator.K0occ", Dimension::P);

    K1    = DS->def<kids_complex>("integrator.K1", Dimension::PFF);
    K1dia = DS->def<kids_complex>("integrator.K1dia", Dimension::PF);
    K1occ = DS->def<kids_complex>("integrator.K1occ", Dimension::P);

    K2    = DS->def<kids_complex>("integrator.K2", Dimension::PFF);
    K2dia = DS->def<kids_complex>("integrator.K2dia", Dimension::PF);
    K2occ = DS->def<kids_complex>("integrator.K2occ", Dimension::P);

    K1QA    = DS->def<kids_complex>("integrator.K1QA", Dimension::PFF);
    K1QAdia = DS->def<kids_complex>("integrator.K1QAdia", Dimension::PF);
    K1QAocc = DS->def<kids_complex>("integrator.K1QAocc", Dimension::P);
    K2QA    = DS->def<kids_complex>("integrator.K2QA", Dimension::PFF);
    K2QAdia = DS->def<kids_complex>("integrator.K2QAdia", Dimension::PF);
    K2QAocc = DS->def<kids_complex>("integrator.K2QAocc", Dimension::P);
    K1DA    = DS->def<kids_complex>("integrator.K1DA", Dimension::PFF);
    K1DAdia = DS->def<kids_complex>("integrator.K1DAdia", Dimension::PF);
    K1DAocc = DS->def<kids_complex>("integrator.K1DAocc", Dimension::P);
    K2DA    = DS->def<kids_complex>("integrator.K2DA", Dimension::PFF);
    K2DAdia = DS->def<kids_complex>("integrator.K2DAdia", Dimension::PF);
    K2DAocc = DS->def<kids_complex>("integrator.K2DAocc", Dimension::P);

    K1QD    = DS->def<kids_complex>("integrator.K1QD", Dimension::PFF);
    K1QDdia = DS->def<kids_complex>("integrator.K1QDdia", Dimension::PF);
    K1QDocc = DS->def<kids_complex>("integrator.K1QDocc", Dimension::P);
    K2QD    = DS->def<kids_complex>("integrator.K2QD", Dimension::PFF);
    K2QDdia = DS->def<kids_complex>("integrator.K2QDdia", Dimension::PF);
    K2QDocc = DS->def<kids_complex>("integrator.K2QDocc", Dimension::P);
    K1DD    = DS->def<kids_complex>("integrator.K1DD", Dimension::PFF);
    K1DDdia = DS->def<kids_complex>("integrator.K1DDdia", Dimension::PF);
    K1DDocc = DS->def<kids_complex>("integrator.K1DDocc", Dimension::P);
    K2DD    = DS->def<kids_complex>("integrator.K2DD", Dimension::PFF);
    K2DDdia = DS->def<kids_complex>("integrator.K2DDdia", Dimension::PF);
    K2DDocc = DS->def<kids_complex>("integrator.K2DDocc", Dimension::P);

    // read input operator
    OpA = DS->def<kids_complex>("integrator.OpA", Dimension::FF);
    OpB = DS->def<kids_complex>("integrator.OpB", Dimension::FF);
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
        kids_complex* K0occ   = this->K0occ + iP;
        kids_complex* K0dia   = this->K0dia + iP * Dimension::F;
        kids_complex* K0      = this->K0 + iP * Dimension::FF;
        kids_complex* K1occ   = this->K1occ + iP;
        kids_complex* K1dia   = this->K1dia + iP * Dimension::F;
        kids_complex* K1      = this->K1 + iP * Dimension::FF;
        kids_complex* K2occ   = this->K2occ + iP;
        kids_complex* K2dia   = this->K2dia + iP * Dimension::F;
        kids_complex* K2      = this->K2 + iP * Dimension::FF;
        kids_complex* K1QAocc = this->K1QAocc + iP;
        kids_complex* K1QAdia = this->K1QAdia + iP * Dimension::F;
        kids_complex* K1QA    = this->K1QA + iP * Dimension::FF;
        kids_complex* K2QAocc = this->K2QAocc + iP;
        kids_complex* K2QAdia = this->K2QAdia + iP * Dimension::F;
        kids_complex* K2QA    = this->K2QA + iP * Dimension::FF;
        kids_complex* K1DAocc = this->K1DAocc + iP;
        kids_complex* K1DAdia = this->K1DAdia + iP * Dimension::F;
        kids_complex* K1DA    = this->K1DA + iP * Dimension::FF;
        kids_complex* K2DAocc = this->K2DAocc + iP;
        kids_complex* K2DAdia = this->K2DAdia + iP * Dimension::F;
        kids_complex* K2DA    = this->K2DA + iP * Dimension::FF;
        kids_complex* K1QDocc = this->K1QDocc + iP;
        kids_complex* K1QDdia = this->K1QDdia + iP * Dimension::F;
        kids_complex* K1QD    = this->K1QD + iP * Dimension::FF;
        kids_complex* K2QDocc = this->K2QDocc + iP;
        kids_complex* K2QDdia = this->K2QDdia + iP * Dimension::F;
        kids_complex* K2QD    = this->K2QD + iP * Dimension::FF;
        kids_complex* K1DDocc = this->K1DDocc + iP;
        kids_complex* K1DDdia = this->K1DDdia + iP * Dimension::F;
        kids_complex* K1DD    = this->K1DD + iP * Dimension::FF;
        kids_complex* K2DDocc = this->K2DDocc + iP;
        kids_complex* K2DDdia = this->K2DDdia + iP * Dimension::F;
        kids_complex* K2DD    = this->K2DD + iP * Dimension::FF;

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
kids_real* Kernel_Elec::T;
kids_real* Kernel_Elec::T_init;

int* Kernel_Elec::occ_nuc;
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
kids_complex* Kernel_Elec::K0occ;
kids_complex* Kernel_Elec::K0dia;
kids_complex* Kernel_Elec::K1;
kids_complex* Kernel_Elec::K1occ;
kids_complex* Kernel_Elec::K1dia;
kids_complex* Kernel_Elec::K2;
kids_complex* Kernel_Elec::K2occ;
kids_complex* Kernel_Elec::K2dia;
kids_complex* Kernel_Elec::K1QA;
kids_complex* Kernel_Elec::K1QAocc;
kids_complex* Kernel_Elec::K1QAdia;
kids_complex* Kernel_Elec::K2QA;
kids_complex* Kernel_Elec::K2QAocc;
kids_complex* Kernel_Elec::K2QAdia;
kids_complex* Kernel_Elec::K1DA;
kids_complex* Kernel_Elec::K1DAocc;
kids_complex* Kernel_Elec::K1DAdia;
kids_complex* Kernel_Elec::K2DA;
kids_complex* Kernel_Elec::K2DAocc;
kids_complex* Kernel_Elec::K2DAdia;
kids_complex* Kernel_Elec::K1QD;
kids_complex* Kernel_Elec::K1QDocc;
kids_complex* Kernel_Elec::K1QDdia;
kids_complex* Kernel_Elec::K2QD;
kids_complex* Kernel_Elec::K2QDocc;
kids_complex* Kernel_Elec::K2QDdia;
kids_complex* Kernel_Elec::K1DD;
kids_complex* Kernel_Elec::K1DDocc;
kids_complex* Kernel_Elec::K1DDdia;
kids_complex* Kernel_Elec::K2DD;
kids_complex* Kernel_Elec::K2DDocc;
kids_complex* Kernel_Elec::K2DDdia;

kids_complex* Kernel_Elec::OpA;
kids_complex* Kernel_Elec::OpB;
kids_complex* Kernel_Elec::TrK1A;
kids_complex* Kernel_Elec::TrK2B;

};  // namespace kids

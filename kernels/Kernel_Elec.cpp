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
    U       = DS->reg<num_complex>("integrator.U", Dimension::PFF);
    c       = DS->reg<num_complex>("integrator.c", Dimension::PF);
    rho_ele = DS->reg<num_complex>("integrator.rho_ele", Dimension::PFF);
    T       = DS->reg<num_real>("model.rep.T", Dimension::PFF);

    occ_nuc = DS->reg<int>("integrator.occ_nuc", Dimension::P);
    rho_nuc = DS->reg<num_complex>("integrator.rho_nuc", Dimension::PFF);

    w = DS->reg<num_complex>("integrator.w");

    K1    = DS->reg<num_complex>("integrator.K1", Dimension::PFF);
    K1dia = DS->reg<num_complex>("integrator.K1dia", Dimension::PF);
    K1occ = DS->reg<num_complex>("integrator.K1occ", Dimension::P);

    K2    = DS->reg<num_complex>("integrator.K2", Dimension::PFF);
    K2dia = DS->reg<num_complex>("integrator.K2dia", Dimension::PF);
    K2occ = DS->reg<num_complex>("integrator.K2occ", Dimension::P);

    K1D    = DS->reg<num_complex>("integrator.K1D", Dimension::PFF);
    K1Ddia = DS->reg<num_complex>("integrator.K1Ddia", Dimension::PF);
    K1Docc = DS->reg<num_complex>("integrator.K1Docc", Dimension::P);

    K2D    = DS->reg<num_complex>("integrator.K2D", Dimension::PFF);
    K2Ddia = DS->reg<num_complex>("integrator.K2Ddia", Dimension::PF);
    K2Docc = DS->reg<num_complex>("integrator.K2Docc", Dimension::P);

    // read input operator
    OpA = DS->reg<num_complex>("integrator.OpA", Dimension::FF);
    OpB = DS->reg<num_complex>("integrator.OpB", Dimension::FF);
}

void Kernel_Elec::init_calc_impl(int stat) {
    exec_kernel_impl(stat);

    double unit = 1.0e0;
    _DataSet->set("init.1", &unit);
    _DataSet->set("init.w", w, Dimension::P);
    _DataSet->set("init.K1", K1, Dimension::PFF);
    _DataSet->set("init.K1dia", K1dia, Dimension::PF);
    _DataSet->set("init.K1occ", K1occ, Dimension::P);
    _DataSet->set("init.K2", K2, Dimension::PFF);
    _DataSet->set("init.K2dia", K2dia, Dimension::PF);
    _DataSet->set("init.K2occ", K2occ, Dimension::P);
    _DataSet->set("init.K1D", K1D, Dimension::PFF);
    _DataSet->set("init.K1Ddia", K1Ddia, Dimension::PF);
    _DataSet->set("init.K1Docc", K1Docc, Dimension::P);
    _DataSet->set("init.K2D", K2D, Dimension::PFF);
    _DataSet->set("init.K2Ddia", K2Ddia, Dimension::PF);
    _DataSet->set("init.K2Docc", K2Docc, Dimension::P);
}

int Kernel_Elec::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* K1occ  = this->K1occ + iP;
        num_complex* K1dia  = this->K1dia + iP * Dimension::F;
        num_complex* K1     = this->K1 + iP * Dimension::FF;
        num_complex* K2occ  = this->K2occ + iP;
        num_complex* K2dia  = this->K2dia + iP * Dimension::F;
        num_complex* K2     = this->K2 + iP * Dimension::FF;
        num_complex* K1Docc = this->K1Docc + iP;
        num_complex* K1Ddia = this->K1Ddia + iP * Dimension::F;
        num_complex* K1D    = this->K1D + iP * Dimension::FF;
        num_complex* K2Docc = this->K2Docc + iP;
        num_complex* K2Ddia = this->K2Ddia + iP * Dimension::F;
        num_complex* K2D    = this->K2D + iP * Dimension::FF;

        ///////////////////////////////////////////////////////

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1dia[i] = K1[ii];
        K1occ[0] = K1[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2dia[i] = K2[ii];
        K2occ[0] = K2[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K1Ddia[i] = K1D[ii];
        K1Docc[0] = K1D[occ0 * Dimension::Fadd1];

        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2Ddia[i] = K2D[ii];
        K2Docc[0] = K2D[occ0 * Dimension::Fadd1];
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
num_real* Kernel_Elec::T;
num_real* Kernel_Elec::T_init;

int* Kernel_Elec::occ_nuc;
num_complex* Kernel_Elec::rho_nuc;
num_complex* Kernel_Elec::rho_nuc_init;

// time correlation function
num_complex* Kernel_Elec::w;
num_complex* Kernel_Elec::K1;
num_complex* Kernel_Elec::K1occ;
num_complex* Kernel_Elec::K1dia;
num_complex* Kernel_Elec::K2;
num_complex* Kernel_Elec::K2occ;
num_complex* Kernel_Elec::K2dia;
num_complex* Kernel_Elec::K1D;
num_complex* Kernel_Elec::K1Docc;
num_complex* Kernel_Elec::K1Ddia;
num_complex* Kernel_Elec::K2D;
num_complex* Kernel_Elec::K2Docc;
num_complex* Kernel_Elec::K2Ddia;

num_complex* Kernel_Elec::OpA;
num_complex* Kernel_Elec::OpB;
num_complex* Kernel_Elec::TrK1A;
num_complex* Kernel_Elec::TrK2B;

};  // namespace PROJECT_NS

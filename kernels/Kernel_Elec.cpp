#include "Kernel_Elec.h"

#include "../core/linalg.h"
#include "Kernel_Dimension.h"
#include "Kernel_Random.h"

namespace PROJECT_NS {

void Kernel_Elec::read_param_impl(Param* PM) {
    occ0 = PM->get<int>("occ", LOC(), -1);
    if (occ0 < 0) throw std::runtime_error("occ < 0");
    if (occ0 >= Kernel_Dimension::F) throw std::runtime_error("occ >= F");
}

void Kernel_Elec::init_data_impl(DataSet* DS) {
    c       = DS->reg<num_complex>("integrator.c", Kernel_Dimension::F);
    rho_ele = DS->reg<num_complex>("integrator.rho_ele", Kernel_Dimension::FF);
    rho_nuc = DS->reg<num_complex>("integrator.rho_nuc", Kernel_Dimension::FF);
    occ_ele = DS->reg<int>("integrator.occ_ele");
    occ_nuc = DS->reg<int>("integrator.occ_nuc");

    w      = DS->reg<num_complex>("integrator.w");
    K0     = DS->reg<num_complex>("integrator.K0", Kernel_Dimension::FF);
    wK0    = DS->reg<num_complex>("integrator.wK0", Kernel_Dimension::FF);
    wK0dia = DS->reg<num_complex>("integrator.wK0dia", Kernel_Dimension::F);
    wK0occ = DS->reg<num_complex>("integrator.wK0occ");

    Kt    = DS->reg<num_complex>("integrator.Kt", Kernel_Dimension::FF);
    Ktdia = DS->reg<num_complex>("integrator.Ktdia", Kernel_Dimension::F);

    mapvar = DS->reg<num_real>("integrator.mapvar", 2 * Kernel_Dimension::F);

    DS->set("init.c", c, Kernel_Dimension::F);
    DS->set("init.rho_ele", rho_ele, Kernel_Dimension::FF);
    DS->set("init.K0", K0, Kernel_Dimension::FF);
    DS->set("init.wK0", wK0, Kernel_Dimension::FF);
    DS->set("init.wK0dia", wK0dia, Kernel_Dimension::F);
    DS->set("init.wK0occ", wK0occ);
}

void Kernel_Elec::init_calc_impl(int stat) {
    for (int i = 0; i < Kernel_Dimension::FF; ++i) wK0[i] = w[0] * K0[i];
    for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) wK0dia[i] = wK0[ii];
    wK0occ[0] = wK0[occ0 * Kernel_Dimension::Fadd1];

    _DataSet->set("init.c", c, Kernel_Dimension::F);
    _DataSet->set("init.rho_ele", rho_ele, Kernel_Dimension::FF);
    _DataSet->set("init.K0", K0, Kernel_Dimension::FF);
    _DataSet->set("init.wK0", wK0, Kernel_Dimension::FF);
    _DataSet->set("init.wK0dia", wK0dia, Kernel_Dimension::F);
    _DataSet->set("init.wK0occ", wK0occ);
}

int Kernel_Elec::exec_kernel_impl(int stat) {
    for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) Ktdia[i] = Kt[ii];
    return 0;
}

/**
 * @brief convert mapping variables to c (electonic amplititude)
 */
int Kernel_Elec::c_from_mapvar(num_complex* c, num_real* mapvar, int fdim) {
    num_real *mapx = mapvar, *mapp = mapvar + fdim;
    for (int i = 0; i < fdim; ++i) { c[i] = phys::math::sqrthalf * (mapx[i] + phys::math::im * mapp[i]); }
    return 0;
}

/**
 * @brief convert c (electonic amplititude) to mapping variables
 */
int Kernel_Elec::mapvar_from_c(num_real* mapvar, num_complex* c, int fdim) {
    num_real *mapx = mapvar, *mapp = mapvar + fdim;
    for (int i = 0; i < fdim; ++i) {
        mapx[i] = phys::math::sqrttwo * std::real(c[i]);
        mapp[i] = phys::math::sqrttwo * std::imag(c[i]);
    }
    return 0;
}

/**
 * @brief convert c (electonic amplititude) to kernel (affine map of the density)
 */
int Kernel_Elec::ker_from_c(num_complex* ker, num_complex* c, num_real xi, num_real gamma, int fdim) {
    ARRAY_OUTER_CONJ2(ker, c, c, fdim, fdim);
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
int Kernel_Elec::ker_from_rho(num_complex* ker, num_complex* rho, num_real xi, num_real gamma, int fdim) {
    for (int i = 0, idx = 0; i < fdim; ++i) {
        for (int j = 0; j < fdim; ++j, ++idx) {
            ker[idx] = xi * rho[idx];
            if (i == j) ker[idx] -= phys::math::iu * gamma;
        }
    }
    return 0;
}

/**
 * @brief sampling mapping variables from focused condition (A & Q as inputs)
 */
int Kernel_Elec::mapvar_focus(num_real* mapvar, int fdim) {
    num_real *mapA = mapvar, *mapQ = mapvar + fdim;  // actions and angles
    num_real *mapx = mapvar, *mapp = mapvar + fdim;  // x and p
    Kernel_Random::rand_uniform(mapQ, fdim, phys::math::twopi);
    for (int i = 0; i < fdim; ++i) {
        // conversion (A,Q) => (x, p)
        mapx[i] = phys::math::sqrttwo * sqrt(mapA[i]) * cos(mapQ[i]);
        mapp[i] = mapx[i] * tan(mapQ[i]);
    }
    return 0;
}

num_complex* Kernel_Elec::w;  // measure of phase point
num_complex* Kernel_Elec::c;
num_real* Kernel_Elec::mapvar;
num_real* Kernel_Elec::mapx;
num_real* Kernel_Elec::mapp;
num_real* Kernel_Elec::mapA;
num_real* Kernel_Elec::mapQ;
int Kernel_Elec::occ0;

// two densities for dynamic
int* Kernel_Elec::occ_ele;
int* Kernel_Elec::occ_nuc;
num_complex* Kernel_Elec::rho_ele;  // electronic density
num_complex* Kernel_Elec::rho_nuc;  // nuclear weighting density

// time correlation function
num_complex* Kernel_Elec::K0;
num_complex* Kernel_Elec::wK0;
num_complex* Kernel_Elec::wK0occ;
num_complex* Kernel_Elec::wK0dia;
num_complex* Kernel_Elec::Kt;
num_complex* Kernel_Elec::Ktdia;

num_complex* Kernel_Elec::OpA;
num_complex* Kernel_Elec::OpB;
num_complex* Kernel_Elec::TrK0A;
num_complex* Kernel_Elec::TrKtB;


};  // namespace PROJECT_NS

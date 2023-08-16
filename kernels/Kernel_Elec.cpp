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
    c       = DS->reg<num_complex>("integrator.c", Dimension::F);
    rho_ele = DS->reg<num_complex>("integrator.rho_ele", Dimension::FF);
    rho_nuc = DS->reg<num_complex>("integrator.rho_nuc", Dimension::FF);
    gmat    = DS->reg<num_complex>("integrator.gmat", Dimension::FF);
    occ_nuc = DS->reg<int>("integrator.occ_nuc");

    w      = DS->reg<num_complex>("integrator.w");
    K1     = DS->reg<num_complex>("integrator.K1", Dimension::FF);
    wK1    = DS->reg<num_complex>("integrator.wK1", Dimension::FF);
    wK1dia = DS->reg<num_complex>("integrator.wK1dia", Dimension::F);
    wK1occ = DS->reg<num_complex>("integrator.wK1occ");

    K2    = DS->reg<num_complex>("integrator.K2", Dimension::FF);
    K2dia = DS->reg<num_complex>("integrator.K2dia", Dimension::F);

    K1Q = DS->reg<num_complex>("integrator.K1Q", Dimension::FF);
    K2Q = DS->reg<num_complex>("integrator.K2Q", Dimension::FF);

    mapvar = DS->reg<num_real>("integrator.mapvar", 2 * Dimension::F);

    double unit = 1.0e0;
    DS->set("init.1", &unit);
    DS->set("init.w", w);
    DS->set("init.c", c, Dimension::F);
    DS->set("init.rho_ele", rho_ele, Dimension::FF);
    DS->set("init.K1", K1, Dimension::FF);
    DS->set("init.wK1", wK1, Dimension::FF);
    DS->set("init.wK1dia", wK1dia, Dimension::F);
    DS->set("init.wK1occ", wK1occ);
    DS->set("init.K1Q", K1Q, Dimension::FF);
    DS->set("init.K2Q", K2Q, Dimension::FF);
    DS->set("init.K2", K2, Dimension::FF);
}

void Kernel_Elec::init_calc_impl(int stat) {
    for (int i = 0; i < Dimension::FF; ++i) wK1[i] = w[0] * K1[i];
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) wK1dia[i] = wK1[ii];
    wK1occ[0] = wK1[occ0 * Dimension::Fadd1];

    _DataSet->set("init.w", w);
    _DataSet->set("init.c", c, Dimension::F);
    _DataSet->set("init.rho_ele", rho_ele, Dimension::FF);
    _DataSet->set("init.K1", K1, Dimension::FF);
    _DataSet->set("init.wK1", wK1, Dimension::FF);
    _DataSet->set("init.wK1dia", wK1dia, Dimension::F);
    _DataSet->set("init.wK1occ", wK1occ);
    _DataSet->set("init.K1Q", K1Q, Dimension::FF);
    _DataSet->set("init.K2Q", K2Q, Dimension::FF);
    _DataSet->set("init.K2", K2, Dimension::FF);
}

int Kernel_Elec::exec_kernel_impl(int stat) {
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) K2dia[i] = K2[ii];
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
 * @brief convert c (electonic amplititude) to kernel (affine map of the density)
 */
int Kernel_Elec::ker_from_rho_quantize(num_complex* ker, num_complex* rho, num_real xi, num_real gamma, int occt,
                                       int fdim) {
    for (int i = 0, idx = 0; i < fdim; ++i) {
        for (int j = 0; j < fdim; ++j, ++idx) {
            ker[idx] = xi * rho[idx];
            if (i == j) ker[idx] = (i == occt) ? 1.0e0 : 0.0e0;
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
int* Kernel_Elec::occ_nuc;
num_complex* Kernel_Elec::rho_ele;  // electronic density
num_complex* Kernel_Elec::rho_nuc;  // nuclear weighting density
num_complex* Kernel_Elec::gmat;     // nuclear weighting density

// time correlation function
num_complex* Kernel_Elec::K1;
num_complex* Kernel_Elec::wK1;
num_complex* Kernel_Elec::wK1occ;
num_complex* Kernel_Elec::wK1dia;
num_complex* Kernel_Elec::K2;
num_complex* Kernel_Elec::K2dia;
num_complex* Kernel_Elec::K1Q;
num_complex* Kernel_Elec::K2Q;

num_complex* Kernel_Elec::OpA;
num_complex* Kernel_Elec::OpB;
num_complex* Kernel_Elec::TrK1A;
num_complex* Kernel_Elec::TrK2B;


};  // namespace PROJECT_NS

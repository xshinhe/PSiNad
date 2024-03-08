#include "Model_ElectronTransfer.h"

// #include <glog/logging.h>

#include "../core/linalg.h"
#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_NADForce.h"
#include "../kernels/Kernel_Random.h"

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

void Model_ElectronTransfer::read_param_impl(Param* PM) {
    // size information
    Nb        = _Param->get<int>("Nb", LOC());
    nbath     = _Param->get<int>("nbath", LOC());
    scan_flag = _Param->get<int>("scan_flag", LOC(), 0);

    // CHECK_EQ(nbath, 1);
    // CHECK_EQ(Nb + 1, Dimension::N);

    FORCE_OPT::nbath = nbath;
    FORCE_OPT::Nb    = Nb;  //@ bugs
    // omega0  = 3.5e-4;
    // lambda0 = 2.39e-2;
    // c0      = std::sqrt(0.5 * lambda0) * omega0;
    // eta = 1.49e-5;
    // alpha = 2.0 * eta / phys::math::pi;
    // omegac = omega0;
    // lambda = 0.5*alpha*omegac
}

void Model_ElectronTransfer::init_data_impl(DataSet* DS) {
    /// 1) System
    Hsys = DS->def<kids_real>("model.Hsys", Dimension::FF);
    memset(Hsys, 0, Dimension::FF * sizeof(kids_real));

    omega0             = _Param->get<double>("omega0", LOC(), 3.5e-4);
    lambda0            = _Param->get<double>("lambda0", LOC(), 2.39e-2);
    coeff0             = _Param->get<double>("coeff0", LOC(), std::sqrt(0.5 * lambda0) * omega0);
    double temperature = _Param->get<double>("temperature", LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman
    double fullbias    = _Param->get<double>("fullbias", LOC(), phys::energy_d, 0.0f);
    double delta       = _Param->get<double>("delta", LOC(), phys::energy_d, 5.0e-5);
    // double eta         = 1.49e-5; // eta has different relation with alpha!!!

    switch (scan_flag) {
        case 1:
            fullbias = 0.0e0;
            beta     = 1.0 / 9.5e-4;
            delta    = delta / beta;
            break;
        case 2:
            fullbias = fullbias * lambda0;
            // beta     = 1.0 / 9.5e-4;
            delta = 5.0e-5;
            break;
        default: {
            break;
        }
    }
    double bias   = fullbias / 2;
    double HSB[4] = {bias, delta, delta, -bias};
    for (int i = 0; i < Dimension::FF; ++i) Hsys[i] = HSB[i];

    /// 2) init Bath sub-kernel (declaration & call)
    for (auto pkernel : _kernel_vector) pkernel->init_data(DS);
    omegas  = DS->def<double>("model.bath.omegas", Nb);
    coeffs  = DS->def<double>("model.bath.coeffs", Nb);
    x_sigma = DS->def<double>("model.bath.x_sigma", Nb);
    p_sigma = DS->def<double>("model.bath.p_sigma", Nb);

    double w2 = omega0 * omega0;
    for (int j = 0; j < Nb; ++j) { w2 += (coeffs[j] * coeffs[j]) / (omegas[j] * omegas[j]); }
    omega0 = std::sqrt(w2);

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Q    = DS->def<double>("model.coupling.Q", nbath * Dimension::FF);
    Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;

    // model field
    mass = DS->def<double>("model.mass", Dimension::N);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;

    vpes = DS->def<double>("model.vpes", Dimension::P);
    grad = DS->def<double>("model.grad", Dimension::PN);
    hess = DS->def<double>("model.hess", Dimension::PNN);
    V    = DS->def<double>("model.V", Dimension::PFF);
    dV   = DS->def<double>("model.dV", Dimension::PNFF);
    // ddV  = DS->def<double>("model.ddV", Dimension::NNFF);

    // init & integrator
    x = DS->def<double>("integrator.x", Dimension::PN);
    p = DS->def<double>("integrator.p", Dimension::PN);

    // ARRAY_SHOW(Hsys, Dimension::F, Dimension::F);
    // ARRAY_SHOW(omegas, 1, Nb);
    // ARRAY_SHOW(coeffs, 1, Nb);
    // ARRAY_SHOW(x_sigma, 1, Nb);
    // ARRAY_SHOW(p_sigma, 1, Nb);
    // exit(-1);
}

void Model_ElectronTransfer::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* x = this->x + iP * Dimension::N;
        kids_real* p = this->p + iP * Dimension::N;

        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);

        x[0] = x[0] * std::sqrt(0.5e0 / (omega0 * tanh(0.5e0 * beta * omega0)))  //
               - coeff0 / (omega0 * omega0);
        p[0] = p[0] * std::sqrt(0.5e0 * omega0 / (tanh(0.5e0 * beta * omega0)));
        for (int j = 0, jadd1 = 1; j < Nb; ++j, ++jadd1) {
            x[jadd1] = x[jadd1] * std::sqrt(0.5e0 / (omegas[j] * tanh(0.5e0 * beta * omegas[j])));
            p[jadd1] = p[jadd1] * std::sqrt(0.5e0 * omegas[j] / (tanh(0.5e0 * beta * omegas[j])));
        }
    }

    _DataSet->def("init.x", x, Dimension::PN);
    _DataSet->def("init.p", p, Dimension::PN);
    exec_kernel(stat);
}

int Model_ElectronTransfer::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* x    = this->x + iP * Dimension::N;
        kids_real* vpes = this->vpes + iP;
        kids_real* grad = this->grad + iP * Dimension::N;
        kids_real* hess = this->hess + iP * Dimension::NN;
        kids_real* V    = this->V + iP * Dimension::FF;
        kids_real* dV   = this->dV + iP * Dimension::NFF;

        // note we implement mass = 1

        // calculate nuclear vpes and grad
        vpes[0] = omega0 * omega0 * x[0] * x[0];
        grad[0] = omega0 * omega0 * x[0];
        for (int j = 0, jadd1 = 1; j < Nb; ++j, ++jadd1) {
            vpes[0] += omegas[j] * omegas[j] * x[jadd1] * x[jadd1];
            grad[jadd1] = omegas[j] * omegas[j] * x[jadd1];
        }
        vpes[0] *= 0.5e0;
        for (int j = 0, jadd1 = 1; j < Nb; ++j, ++jadd1) {
            vpes[0] += coeffs[j] * x[0] * x[jadd1];
            grad[0] += coeffs[j] * x[jadd1];
            grad[jadd1] += coeffs[j] * x[0];
        }

        // calculate electronic V and dV
        for (int i = 0; i < Dimension::FF; ++i) V[i] = Hsys[i];
        V[0] += coeff0 * x[0];
        V[3] -= coeff0 * x[0];

        if (count_exec == 0) {  // only calculate once time
            ARRAY_CLEAR(dV, Dimension::NFF);
            dV[0] = coeff0;
            dV[3] = -coeff0;
        }
    }
    return 0;
}

};  // namespace PROJECT_NS

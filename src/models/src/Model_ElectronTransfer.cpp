#include "kids/Model_ElectronTransfer.h"

#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_ElectronTransfer::getName() { return "Model_ElectronTransfer"; }

int Model_ElectronTransfer::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_ElectronTransfer::setInputParam_impl(std::shared_ptr<Param> PM) {
    // size information
    Nb        = _param->get_int({"model.Nb"}, LOC());
    nbath     = _param->get_int({"model.nbath"}, LOC());
    scan_flag = _param->get_int({"model.scan_flag"}, LOC(), 0);

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

void Model_ElectronTransfer::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    /// 1) System
    Hsys = DS->def_real("model.Hsys", Dimension::FF);
    memset(Hsys, 0, Dimension::FF * sizeof(kids_real));

    omega0             = _param->get_real({"model.omega0"}, LOC(), 3.5e-4);
    lambda0            = _param->get_real({"model.lambda0"}, LOC(), 2.39e-2);
    coeff0             = _param->get_real({"model.coeff0"}, LOC(), std::sqrt(0.5 * lambda0) * omega0);
    double temperature = _param->get_real({"model.temperature"}, LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman
    double fullbias    = _param->get_real({"model.fullbias"}, LOC(), phys::energy_d, 0.0f);
    double delta       = _param->get_real({"model.delta"}, LOC(), phys::energy_d, 5.0e-5);
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
    for (auto pkernel : _child_kernels) pkernel->setInputDataSet(DS);
    omegas  = DS->def(DATA::model::bath::omegas);
    coeffs  = DS->def(DATA::model::bath::coeffs);
    x_sigma = DS->def(DATA::model::bath::x_sigma);
    p_sigma = DS->def(DATA::model::bath::p_sigma);

    double w2 = omega0 * omega0;
    for (int j = 0; j < Nb; ++j) { w2 += (coeffs[j] * coeffs[j]) / (omegas[j] * omegas[j]); }
    omega0 = std::sqrt(w2);

    /// 3) bilinear Coupling (saving order: L, nbath, Nb, FF)
    Q    = DS->def(DATA::model::coupling::Q);
    Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;

    // model field
    mass = DS->def(DATA::model::mass);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;

    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    hess = DS->def(DATA::model::hess);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    // ddV  = DS->def(DATA::model::ddV);

    // init & integrator
    x = DS->def(DATA::integrator::x);
    p = DS->def(DATA::integrator::p);
}

Status& Model_ElectronTransfer::initializeKernel_impl(Status& stat) {
    executeKernel(stat);
    return stat;  // @todo

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

    _dataset->def_real("init.x", x, Dimension::PN);
    _dataset->def_real("init.p", p, Dimension::PN);
    executeKernel(stat);
    return stat;
}

Status& Model_ElectronTransfer::executeKernel_impl(Status& stat) {
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
    return stat;
}

};  // namespace PROJECT_NS

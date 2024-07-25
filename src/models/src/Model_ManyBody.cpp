#include "kids/Model_ManyBody.h"

#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_ManyBody::getName() { return "Model_ManyBody"; }

int Model_ManyBody::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_ManyBody::setInputParam_impl(std::shared_ptr<Param> PM) {
    manybody_type = ManyBodyPolicy::_from(_param->get_string({"model.manybody_flag"}, LOC(), "Ising"));
    Jp            = PM->get_real({"model.Jp"}, LOC(), 0.0e0);
    Jz            = PM->get_real({"model.Jz"}, LOC(), 1.0e0);
    alpha         = PM->get_real({"model.alpha"}, LOC(), 0.0e0);
    omega         = PM->get_real({"model.omega"}, LOC(), 0.0e0);
}

void Model_ManyBody::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    mass = DS->def(DATA::model::mass);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;
    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    hess = DS->def(DATA::model::hess);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    // ddV  = DS->def(DATA::model::ddV);
    w       = DS->def(DATA::model::w);
    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);
    x       = DS->def(DATA::integrator::x);
    p       = DS->def(DATA::integrator::p);

    // for many body interaction
    Jpmat   = DS->def(DATA::model::MB::Jpmat);
    Jzmat   = DS->def(DATA::model::MB::Jzmat);
    SXred   = DS->def(DATA::model::MB::SXred);
    SYred   = DS->def(DATA::model::MB::SYred);
    SZred   = DS->def(DATA::model::MB::SZred);
    H1      = DS->def(DATA::model::MB::H1);
    H2      = DS->def(DATA::model::MB::H2);
    H       = DS->def(DATA::model::rep::H);
    rho_ele = DS->def(DATA::integrator::rho_ele);

    /** Ising model and XY model
        Jij ≡ J*[1 − 3cos^2(θ)] / rij^α
    */
    for (int i = 0, ik = 0; i < Dimension::P; ++i) {
        for (int k = 0; k < Dimension::P; ++k, ++ik) {
            if (i == k) {
                Jpmat[ik] = 0.0f, Jzmat[ik] = 0.0f;
                continue;
            }
            int dist  = i > k ? i - k : k - i;  // or not? @data-202407
            Jpmat[ik] = Jp / std::pow(dist, alpha);
            Jzmat[k]  = Jz / std::pow(dist, alpha);
        }
    }
}

kids_complex PauliX[4] = {phys::math::iz, phys::math::iu, phys::math::iu, phys::math::iz};
kids_complex PauliY[4] = {phys::math::iz, -phys::math::im, phys::math::im, phys::math::iz};
kids_complex PauliZ[4] = {phys::math::iu, phys::math::iz, phys::math::iz, -phys::math::iu};

Status& Model_ManyBody::execute_Heisenberg(Status& stat) {
    memset(H, 0, Dimension::PFF * sizeof(kids_complex));
    // for (int i = 0; i < Dimension::PFF; ++i) H[i] = H0[i];  // constant Hamiltonian
    // reduce X/Y/Z in temporary array redX, redY, redZ (saving time)
    for (int i = 0; i < Dimension::P; ++i) {
        kids_complex* rhoi = rho_ele + i * Dimension::FF;
        SXred[i]           = ARRAY_TRACE2(rhoi, PauliX, Dimension::F, Dimension::F);
        SYred[i]           = ARRAY_TRACE2(rhoi, PauliY, Dimension::F, Dimension::F);
        SZred[i]           = ARRAY_TRACE2(rhoi, PauliZ, Dimension::F, Dimension::F);
    }
    for (int i = 0, ij = 0; i < Dimension::P; ++i) {
        kids_complex* Hi = H + i * Dimension::FF;

        // one-body interaction
        for (int k = 0; k < Dimension::FF; ++k) Hi[k] = omega * PauliX[k];

        // two-body interaction
        for (int j = 0; j < Dimension::P; ++j, ++ij) {
            for (int k = 0; k < Dimension::FF; ++k) Hi[k] += 0.25e0 * PauliX[k] * Jpmat[ij] * SXred[j];
            for (int k = 0; k < Dimension::FF; ++k) Hi[k] += 0.25e0 * PauliY[k] * Jpmat[ij] * SYred[j];
            for (int k = 0; k < Dimension::FF; ++k) Hi[k] += 0.5e0 * PauliZ[k] * Jzmat[ij] * SZred[j];
            // three-body interaction: pass
            // ...
        }
    }
    return stat;
};

Status& Model_ManyBody::executeKernel_impl(Status& stat) {
    switch (manybody_type) {
        case ManyBodyPolicy::TwoSite: {
            break;
        }
        case ManyBodyPolicy::Ising:
        case ManyBodyPolicy::XY:
        case ManyBodyPolicy::Heisenberg: {
            execute_Heisenberg(stat);
            break;
        }
        case ManyBodyPolicy::Dicke:
        case ManyBodyPolicy::TavisCummings: {
            break;
        }
        case ManyBodyPolicy::Anderson:
        case ManyBodyPolicy::Hubbard: {
            break;
        }
        case ManyBodyPolicy::BowTie: {
            break;
        }
    }
    return stat;
}

};  // namespace PROJECT_NS

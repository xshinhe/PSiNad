#include <cassert>
#include <iostream>

#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Representation.h"
#include "kids/Status.h"
#include "kids/vars_list.h"

using namespace PROJECT_NS;

static void init_minimal_dims(std::size_t M = 1, std::size_t P = 1, std::size_t N = 1, std::size_t F = 2) {
    Dimension::M = M;
    Dimension::P = P;
    Dimension::N = N;
    Dimension::F = F;
    Dimension::static_build_shapes();
}

int main() {
    // Arrange: minimal 2-level system, 1 nuclear DOF
    init_minimal_dims(1, 1, 1, 2);

    // Keep both input and nuclear representations in Diabatic to avoid extra transforms
    Kernel_Representation::inp_repr_type = RepresentationPolicy::Diabatic;
    Kernel_Representation::nuc_repr_type = RepresentationPolicy::Diabatic;

    // Minimal parameter set
    Param PM("{\n"
             "  \"solver.naforce\": \"EHR\",\n"
             "  \"solver.BATH_FORCE_BILINEAR\": false,\n"
             "  \"solver.offd_projected\": true\n"
             "}",
             Param::fromString);

    // DataSet and variables
    auto DS = std::make_shared<DataSet>();

    // Control
    auto dt = DS->def(DATA::control::dt); dt[0] = 0.1;
    auto succ = DS->def(DATA::control::succ); succ[0] = 1;

    // Nuclear arrays
    auto f    = DS->def(DATA::integrator::f);      for (std::size_t j = 0; j < Dimension::N; ++j) f[j] = 0.0;
    auto fadd = DS->def(DATA::integrator::fadd);   for (std::size_t j = 0; j < Dimension::N; ++j) fadd[j] = 0.0;
    auto p    = DS->def(DATA::integrator::p);      for (std::size_t j = 0; j < Dimension::N; ++j) p[j] = 0.0;
    auto m    = DS->def(DATA::integrator::m);      for (std::size_t j = 0; j < Dimension::N; ++j) m[j] = 1.0;
    auto grad = DS->def(DATA::model::grad);        for (std::size_t j = 0; j < Dimension::N; ++j) grad[j] = 0.0;

    // Electronic Hamiltonian (diabatic) and derivative
    auto V  = DS->def(DATA::model::V);   // 2x2
    auto dV = DS->def(DATA::model::dV);  // N x 2x2, here N=1 so just one 2x2 block

    // Set V = [[0, 0],[0, 1]]
    for (std::size_t ik = 0; ik < Dimension::FF; ++ik) V[ik] = 0.0;
    V[0 * Dimension::F + 0] = 0.0; // (0,0)
    V[1 * Dimension::F + 1] = 1.0; // (1,1)

    // Set dV (for j=0): only (0,0) = 1.5, others 0
    for (std::size_t ik = 0; ik < Dimension::FF; ++ik) dV[ik] = 0.0;
    dV[0 * Dimension::F + 0] = 1.5; // dV00

    // Representation matrices
    auto T = DS->def(DATA::model::rep::T); // identity 2x2
    for (std::size_t ik = 0; ik < Dimension::FF; ++ik) T[ik] = 0.0;
    T[0 * Dimension::F + 0] = 1.0;
    T[1 * Dimension::F + 1] = 1.0;

    // Rho (nuclear and electronic)
    auto rho_ele = DS->def(DATA::integrator::rho_ele);
    auto rho_nuc = DS->def(DATA::integrator::rho_nuc);
    for (std::size_t ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = kids_complex(0.0, 0.0);
    for (std::size_t ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = kids_complex(0.0, 0.0);
    // rho_nuc = |0><0|
    rho_nuc[0 * Dimension::F + 0] = kids_complex(1.0, 0.0);

    // Occupation (not used in EHR branch for force formula, but defined for completeness)
    auto occ_nuc = DS->def(DATA::integrator::occ_nuc);
    occ_nuc[0] = 0;

    // Other required buffers
    auto dE   = DS->def(DATA::model::rep::dE);     for (std::size_t i = 0; i < Dimension::NFF; ++i) dE[i] = 0.0;
    auto vpes = DS->def(DATA::model::vpes);        vpes[0] = 0.0;
    auto Epot = DS->def(DATA::integrator::Epot);   Epot[0] = 0.0;
    auto alpha= DS->def(DATA::integrator::alpha);  alpha[0] = 0.0;
    auto ftmp = DS->def(DATA::integrator::tmp::ftmp);
    auto fproj= DS->def(DATA::integrator::tmp::fproj);
    auto wrho = DS->def(DATA::integrator::tmp::wrho);

    // Act: run the kernel
    Kernel_NAForce knaf;
    knaf.setInputParam(std::make_shared<Param>(PM));
    knaf.setInputDataSet(DS);

    Status stat;
    knaf.initializeKernel(stat); // init calls execute once internally

    // Assert
    if (std::abs(f[0] - 1.5) > 1e-12) {
        std::cerr << "Test failed: f[0]=" << f[0] << ", expected 1.5\n";
        return 1;
    }
    std::cout << "Test passed: f[0]=" << f[0] << "\n";
    return 0;
}

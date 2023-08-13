#include "Kernel_Elec_SH.h"

#include "Kernel_Dimension.h"
#include "Kernel_NADForce.h"
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

int Kernel_Elec_SH::max_choose(num_complex* rho) {
    int imax      = 0;
    num_real vmax = 0.0f;
    for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) {
        if (std::abs(rho[ii]) > vmax) {
            vmax = std::abs(rho[ii]);
            imax = i;
        }
    }
    return imax;
}

int Kernel_Elec_SH::pop_choose(num_complex* rho) {
    num_real rand_tmp;
    num_real sum = 0.0f;
    Kernel_Random::rand_uniform(&rand_tmp);
    for (int i = 0, ii = 0; i < Kernel_Dimension::F; ++i, ii += Kernel_Dimension::Fadd1) {
        sum += std::abs(rho[ii]);
        if (rand_tmp < sum) return i;
    }
    return 0;
}

int Kernel_Elec_SH::hopping_choose(num_complex* rho, num_complex* H, int from, num_real dt) {
    int to = from;
    num_real rand_tmp, sumprob = 0.0f;
    num_real rhoii = std::real(rho[from * Kernel_Dimension::Fadd1]);
    Kernel_Random::rand_uniform(&rand_tmp);

    for (int n = 0; n < Kernel_Dimension::F; ++n) {
        num_real prob =
            (n == from) ? 0.0f
                        : -2.0f * std::imag(rho[n * Kernel_Dimension::F + from] * H[from * Kernel_Dimension::F + n]) /
                              rhoii * dt;
        prob = (prob > 1.0f) ? 1.0f : ((prob < 0.0f) ? 0.0f : prob);  // hopping cut-off
        sumprob += prob;
        if (rand_tmp < sumprob) {
            to = n;
            break;
        }
    }
    return to;
}

void Kernel_Elec_SH::hopping_direction(num_real* direction, num_real* dE, int from, int to) {
    if (to == from) return;
    for (int i = 0; i < Kernel_Dimension::N; ++i) {
        direction[i] = dE[i * Kernel_Dimension::FF + from * Kernel_Dimension::F + to];
    }
}

int Kernel_Elec_SH::hopping_impulse(num_real* direction, num_real* np, num_real* nm, num_real* E,  //
                                    int from, int to, bool reflect) {
    if (to == from) return from;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    num_real coeffa = 0.0f, coeffb = 0.0f, coeffc = E[to] - E[from];
    for (int i = 0; i < Kernel_Dimension::N; ++i) {
        coeffa += 0.5f * direction[i] * direction[i] / nm[i];
        coeffb += np[i] / nm[i] * direction[i];
    }
    coeffb /= coeffa, coeffc /= coeffa;  // normalization for safety

    num_real coeffd = coeffb * coeffb - 4 * coeffc;
    if (coeffd > 0) {
        num_real x1 = 0.5f * (-coeffb + sqrt(coeffd)), x2 = 0.5f * (-coeffb - sqrt(coeffd));
        num_real xx = (std::abs(x1) < std::abs(x2)) ? x1 : x2;
        for (int i = 0; i < Kernel_Dimension::N; ++i) np[i] += xx * direction[i];
        return to;
    } else if (reflect) {  // 2008Algorithm
        num_real xx = -coeffb;
        for (int i = 0; i < Kernel_Dimension::N; ++i) np[i] += xx * direction[i];
        return from;
    } else {  // 1990Algorithm, do nothing
        return from;
    }
    return from;
}


void Kernel_Elec_SH::read_param_impl(Param* PM) {
    sh_type = SHPolicy::_from(PM->get<std::string>("sh_flag", LOC(), "FSSH"));
    reflect = PM->get<bool>("reflect", LOC(), false);
    dt      = PM->get<double>("dt", LOC(), phys::time_d);

    std::string REP_CHECK = PM->get<std::string>("representation_flag", LOC(), "Diabatic");
    assert(REP_CHECK == "Adiabatic");
}

void Kernel_Elec_SH::init_data_impl(DataSet* DS) {
    x         = DS->reg<num_real>("integrator.x", Kernel_Dimension::N);
    p         = DS->reg<num_real>("integrator.p", Kernel_Dimension::N);
    m         = DS->reg<num_real>("integrator.m", Kernel_Dimension::N);
    direction = DS->reg<num_real>("integrator.direction", Kernel_Dimension::N);
    E         = DS->reg<num_real>("model.rep.E", Kernel_Dimension::F);
    dE        = DS->reg<num_real>("model.rep.dE", Kernel_Dimension::NFF);
    T         = DS->reg<num_real>("model.rep.T", Kernel_Dimension::FF);
    H         = DS->reg<num_complex>("model.rep.H", Kernel_Dimension::FF);
}

void Kernel_Elec_SH::init_calc_impl(int stat) {
    Kernel_NADForce::NADForce_type = NADForcePolicy::BO;  // set BO dynamics for SH

    Kernel_Elec::w[0] = 1.0e0;

    for (int i = 0, ik = 0; i < Kernel_Dimension::F; ++i) {
        for (int k = 0; k < Kernel_Dimension::F; ++k, ++ik) {
            Kernel_Elec::rho_ele[ik] = (i == Kernel_Elec::occ0 && k == Kernel_Elec::occ0) ? 1.0e0 : 0.0e0;
        }
    }

    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Kernel_Dimension::F, Kernel_Dimension::F,
                             Kernel_Dimension::F, Kernel_Dimension::F);
    }
    *Kernel_Elec::occ_nuc = pop_choose(Kernel_Elec::rho_ele);
    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Kernel_Dimension::F, Kernel_Dimension::F,
                             Kernel_Dimension::F, Kernel_Dimension::F);
    }

    for (int ik = 0; ik < Kernel_Dimension::FF; ++ik) Kernel_Elec::rho_nuc[ik] = Kernel_Elec::rho_ele[ik];
    for (int ik = 0; ik < Kernel_Dimension::FF; ++ik) Kernel_Elec::K0[ik] = Kernel_Elec::rho_ele[ik];
    exec_kernel(stat);
}

int Kernel_Elec_SH::exec_kernel_impl(int stat) {
    /**
     * additional evolution appended to Kernel_Update_rho (auxilrary variables)
     */

    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Kernel_Dimension::F, Kernel_Dimension::F,
                             Kernel_Dimension::F, Kernel_Dimension::F);
    }

    // step 1: determine where to hop
    int to = hopping_choose(Kernel_Elec::rho_ele, H, *Kernel_Elec::occ_nuc, dt);
    // step 2: determine direction to hop
    hopping_direction(direction, dE, *Kernel_Elec::occ_nuc, to);
    // step 3: try hop
    *Kernel_Elec::occ_nuc = hopping_impulse(direction, p, m, E, *Kernel_Elec::occ_nuc, to, reflect);

    Kernel_Elec::ker_from_rho_quantize(Kernel_Elec::Kt, Kernel_Elec::rho_ele, 1, 0, *Kernel_Elec::occ_nuc,
                                       Kernel_Dimension::F);

    if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T, Kernel_Dimension::F, Kernel_Dimension::F,
                             Kernel_Dimension::F, Kernel_Dimension::F);
        ARRAY_MATMUL3_TRANS2(Kernel_Elec::Kt, T, Kernel_Elec::Kt, T, Kernel_Dimension::F, Kernel_Dimension::F,
                             Kernel_Dimension::F, Kernel_Dimension::F);
    }
    return stat;
}


};  // namespace PROJECT_NS

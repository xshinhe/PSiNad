#include "Kernel_Hopping.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
#include "Kernel_NADForce.h"
#include "Kernel_Random.h"
#include "Kernel_Representation.h"

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

int Kernel_Hopping::max_choose(kids_complex* rho) {
    int imax       = 0;
    kids_real vmax = 0.0f;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        if (std::real(rho[ii]) > vmax) {
            vmax = std::real(rho[ii]);
            imax = i;
        }
    }
    return imax;
}

int Kernel_Hopping::pop_choose(kids_complex* rho) {
    kids_real rand_tmp;
    kids_real sum = 0.0f;
    Kernel_Random::rand_uniform(&rand_tmp);
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        sum += std::min({std::max({std::real(rho[ii]), 0.0e0}), 1.0e0});
        if (rand_tmp < sum) return i;
    }
    return 0;
}

int Kernel_Hopping::pop_neg_choose(kids_complex* rho) {
    kids_real rand_tmp;
    kids_real sum = 0.0f;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) sum += std::abs(rho[ii]);
    Kernel_Random::rand_uniform(&rand_tmp);
    rand_tmp *= sum;
    sum = 0.0e0;
    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
        sum += std::abs(rho[ii]);
        if (rand_tmp < sum) return i;
    }
    return 0;
}

int Kernel_Hopping::hopping_choose(kids_complex* rho, kids_complex* H, int from, kids_real dt) {
    int to = from;
    kids_real rand_tmp, sumprob = 0.0f;
    kids_real rhoii = std::real(rho[from * Dimension::Fadd1]);
    Kernel_Random::rand_uniform(&rand_tmp);

    for (int n = 0; n < Dimension::F; ++n) {
        kids_real prob =
            (n == from) ? 0.0f
                        : -2.0f * std::imag(rho[n * Dimension::F + from] * H[from * Dimension::F + n]) / rhoii * dt;
        prob = (prob > 1.0f) ? 1.0f : ((prob < 0.0f) ? 0.0f : prob);  // hopping cut-off
        sumprob += prob;
        if (rand_tmp < sumprob) {
            to = n;
            break;
        }
    }
    return to;
}

void Kernel_Hopping::hopping_direction(kids_real* direction, kids_real* dE, int from, int to) {
    if (to == from) return;
    for (int i = 0; i < Dimension::N; ++i) { direction[i] = dE[i * Dimension::FF + from * Dimension::F + to]; }
}

int Kernel_Hopping::hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm, kids_real* E,  //
                                    int from, int to, bool reflect) {
    if (to == from) return from;

    // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
    kids_real coeffa = 0.0f, coeffb = 0.0f, coeffc = E[to] - E[from];
    for (int i = 0; i < Dimension::N; ++i) {
        coeffa += 0.5f * direction[i] * direction[i] / nm[i];
        coeffb += np[i] / nm[i] * direction[i];
    }
    coeffb /= coeffa, coeffc /= coeffa;  // normalization for safety

    kids_real coeffd = coeffb * coeffb - 4 * coeffc;
    if (coeffd > 0) {
        kids_real x1 = 0.5f * (-coeffb + sqrt(coeffd)), x2 = 0.5f * (-coeffb - sqrt(coeffd));
        kids_real xx = (std::abs(x1) < std::abs(x2)) ? x1 : x2;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return to;
    } else if (reflect) {  // 2008Algorithm
        kids_real xx = -coeffb;
        for (int i = 0; i < Dimension::N; ++i) np[i] += xx * direction[i];
        return from;
    } else {  // 1990Algorithm, do nothing
        return from;
    }
    return from;
}

};  // namespace PROJECT_NS

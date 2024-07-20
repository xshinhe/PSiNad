/**
 * @file Kernel_Elec.h
 * @author xshinhe
 * @version 1.1
 * @date 2023-03
 * @brief initialization kernels for electonic DOFs
 * @details
 *  The initialization of electonic DOFs are tightly related to Solver.
 *  Use it in Solver's Kernel_Builder();
 */

#ifndef Kernel_Elec_Utils_H
#define Kernel_Elec_Utils_H

#include "kids/Kernel.h"
#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/Sampling_Elec.h"
#include "kids/linalg.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

DEFINE_POLICY(SQCPolicy,
              SQR,  // square window
              TRI,  // triangle window
              SPX,  // simplex window
              BIG   // simplex window
);

class elec_utils {
   public:
    /// @{
    static inline double gamma_wigner(int fdim) { return (sqrt(fdim + 1.0e0) - 1) / fdim; }

    static inline double gamma_opt(int fdim) {
        double sum = 0.0e0;
        for (int i = 0; i < fdim; ++i) sum += 1.0e0 / (i + 1);
        return ((fdim - sum) / (sum - 1.0e0)) / fdim;
    }
    /// @}

    /**
     * @brief convert c (electonic amplititude) to kernel (affine map of the density)
     */
    static int ker_from_c(kids_complex* ker, kids_complex* c, kids_real xi, kids_real gamma, int fdim) {
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
    static int ker_from_rho(kids_complex* ker, kids_complex* rho, kids_real xi, kids_real gamma, int fdim,
                            bool quantize = false, int occ = -1) {
        for (int i = 0, ij = 0; i < fdim; ++i) {
            for (int j = 0; j < fdim; ++j, ++ij) {
                ker[ij] = xi * rho[ij];
                if (i == j) ker[ij] = (quantize) ? (i == occ ? 1.0e0 : 0.0e0) : (ker[ij] - gamma);
            }
        }
        return 0;
    }

    static int max_choose(kids_complex* rho) {
        int       imax = 0;
        kids_real vmax = 0.0f;
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            if (std::real(rho[ii]) > vmax) {
                vmax = std::real(rho[ii]);
                imax = i;
            }
        }
        return imax;
    }

    static int pop_choose(kids_complex* rho) {
        kids_real rand_tmp;
        kids_real sum = 0.0f;
        Kernel_Random::rand_uniform(&rand_tmp);
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            sum += std::min({std::max({std::real(rho[ii]), 0.0e0}), 1.0e0});
            if (rand_tmp < sum) return i;
        }
        return 0;
    }

    static int pop_neg_choose(kids_complex* rho) {
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

    static int hopping_choose(kids_complex* rho, kids_complex* H, int from, kids_real dt) {
        int       to = from;
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

    static double calc_ElectricalEnergy(kids_real* E, kids_complex* wrho, int occ) {
        double Ecalc = 0.0e0;
        switch (Kernel_NAForce::NAForce_type) {
            case NAForcePolicy::BO:
            case NAForcePolicy::NAF: {
                Ecalc = E[occ * Dimension::Fadd1];
                break;
            }
            default: {  // EHR, MIX, SD (Eto == Efrom will skip hopping procedure)
                Ecalc = std::real(ARRAY_TRACE2(wrho, E, Dimension::F, Dimension::F));
                break;
            }
        }
        return Ecalc;
    }

    static int calc_distorted_rho(kids_complex* wrho,   // distorted density
                                  kids_complex* rho,    // rho_ele
                                  double        xi,     // xi must be 1
                                  double        gamma,  // gamma must be 0
                                  double        alpha) {
        // initialize distorted-density
        kids_real    L    = 1.0e0 - log(std::abs(alpha));  // @NOTE
        kids_complex norm = 0.0e0;
        for (int i = 0, ik = 0; i < Dimension::F; ++i) {
            for (int k = 0; k < Dimension::F; ++k, ++ik) {
                if (i == k) {
                    wrho[ik] = copysign(1.0, std::real(rho[ik])) * std::pow(std::abs(rho[ik]), L);
                    norm += wrho[ik];
                } else {
                    wrho[ik] = xi * rho[ik];
                }
            }
        }
        // normalization of distorted-density
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) wrho[ii] /= norm;
        return 0;
    }

    static int calc_distorted_force(kids_real*    f1,    // to be calculated
                                    kids_real*    E,     // (input)
                                    kids_real*    dE,    // (input)
                                    kids_complex* wrho,  // distorted rho
                                    kids_complex* rho,   // rho_ele
                                    double        alpha) {
        // initialize distorted-density
        kids_real L            = 1.0e0 - log(std::abs(alpha));  // @NOTE
        kids_real rate_default = (L == 1.0e0) ? 1.0e0 : 0.0e0;

        // averaged adiabatic energy by distorted-density
        double Ew = std::real(ARRAY_TRACE2(wrho, E, Dimension::F, Dimension::F));
        // distortion force (when L= 1 [EHR], the distortion force is zero)
        for (int j = 0, jFF = 0; j < Dimension::N; ++j, jFF += Dimension::FF) {
            kids_real* dEj = dE + jFF;
            f1[j]          = 0.0e0;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                double rate  = ((std::real(rho[ii]) == 0.0e0) ? rate_default : std::real(wrho[ii] / rho[ii]));
                double coeff = (E[ii] - (E[ii] - Ew) * L * rate);
                for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1) {
                    if (i == k) continue;
                    f1[j] +=
                        coeff * std::real(dEj[i * Dimension::F + k] / (E[kk] - E[ii]) * wrho[k * Dimension::F + i] -
                                          dEj[k * Dimension::F + i] / (E[ii] - E[kk]) * wrho[i * Dimension::F + k]);
                }
            }
        }
        return 0;
    }

    static void hopping_direction(kids_real* direction, kids_real* dE, int from, int to) {
        if (to == from) return;
        for (int i = 0; i < Dimension::N; ++i) { direction[i] = dE[i * Dimension::FF + from * Dimension::F + to]; }
    }

    static void hopping_direction(kids_real* direction, kids_real* E, kids_real* dE, kids_complex* rho, int from,
                                  int to) {
        if (to == from) return;

        for (int i = 0; i < Dimension::N; ++i) {
            direction[i] = 0.0f;
            for (int k = 0; k < Dimension::F; ++k)
                if (k != from)
                    direction[i] +=
                        std::real(rho[from * Dimension::F + k] * dE[i * Dimension::FF + k * Dimension::F + from]) /
                        (E[from * Dimension::Fadd1] - E[k * Dimension::Fadd1]);
            for (int k = 0; k < Dimension::F; ++k)
                if (k != to)
                    direction[i] -=
                        std::real(rho[to * Dimension::F + k] * dE[i * Dimension::FF + k * Dimension::F + to]) /
                        (E[to * Dimension::Fadd1] - E[k * Dimension::Fadd1]);
        }
    }

    // static double calc_ElectricalEnergy(kids_real* E, kids_complex* wrho, int occ) {
    //     double Ecalc = 0.0e0;
    //     switch (Kernel_NADForce::NADForce_type) {
    //         case NADForcePolicy::BO:
    //         case NADForcePolicy::CV: {
    //             Ecalc = E[occ * Dimension::Fadd1];
    //             break;
    //         }
    //         default: {  // EHR, MIX, SD (Eto == Efrom will skip hopping procedure)
    //             Ecalc = ARRAY_TRACE2(wrho, E, Dimension::F, Dimension::F);
    //             break;
    //         }
    //     }
    //     return Ecalc;
    // }

    static int hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm,  //
                               kids_real Efrom, kids_real Eto, int from, int to, bool reflect) {
        if (to == from) return from;

        // solve x: Ef + P**2 / (2*M) = Et + (P + direction*x)**2 / (2*M)
        kids_real coeffa = 0.0f, coeffb = 0.0f, coeffc = Eto - Efrom;
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

    static int hopping_impulse(kids_real* direction, kids_real* np, kids_real* nm, kids_real* eig,  //
                               int from, int to, bool reflect) {
        return hopping_impulse(direction, np, nm, eig[from], eig[to], from, to, reflect);
    }

    /**
     * @brief sampling mapping variables from uniform sphere distribution (i.e. uniform simplex for action)
     */
    static int c_sphere(kids_complex* c, int fdim) {
        Kernel_Random::rand_gaussian(reinterpret_cast<double*>(c), 2 * fdim);
        double l = std::sqrt(std::real(ARRAY_INNER_TRANS1(c, c, fdim)));
        for (int i = 0; i < fdim; ++i) c[i] /= l;
        return 0;
    }

    static int c_focus(kids_complex* c, double xi, double gamma, int occ, int fdim) {
        for (int i = 0; i < fdim; ++i) c[i] = ((i == occ ? 1.0e0 : 0.0e0) + gamma) / xi;
        double randu;
        for (int i = 0; i < fdim; ++i) {
            Kernel_Random::rand_uniform(&randu);
            randu *= phys::math::twopi;
            c[i] = std::sqrt(std::abs(c[i])) * (cos(randu) + phys::math::im * sin(randu));
        }
        return 0;
    };

    static int rho_focus(kids_complex* rho, int iocc, double gamma_ou, double gamma_uu, int fdim, bool rand_act,
                         bool pure_phase, bool cont_phase) {
        double randu = 1.0e0;
        if (pure_phase) {
            for (int j = 0; j < fdim; ++j) {
                if (j == iocc) {
                    rho[iocc * fdim + j] = 1.0e0;
                    continue;
                }
                Kernel_Random::rand_uniform(&randu);
                randu = (cont_phase) ? phys::math::twopi * randu : phys::math::halfpi * (int(randu / 0.25f) + 1);
                rho[iocc * fdim + j] = cos(randu) + phys::math::im * sin(randu);
                rho[j * fdim + iocc] = std::conj(rho[iocc * fdim + j]);
            }
            for (int i = 0, ij = 0; i < fdim; ++i) {
                for (int j = 0; j < fdim; ++j, ++ij) {
                    if (i == iocc || j == iocc) continue;
                    rho[ij] = rho[iocc * fdim + j] / rho[iocc * fdim + i];
                }
            }
        } else {
            for (int i = 0; i < fdim; ++i) {
                for (int j = i + 1; j < fdim; ++j) {
                    Kernel_Random::rand_uniform(&randu);
                    randu = (cont_phase) ? phys::math::twopi * randu : phys::math::halfpi * (int(randu / 0.25f) + 1);
                    rho[i * fdim + j] = cos(randu) + phys::math::im * sin(randu);
                    rho[j * fdim + i] = std::conj(rho[i * fdim + j]);
                }
            }
        }
        Kernel_Random::rand_uniform(&randu);
        double occrand = (rand_act) ? int(randu * fdim) : iocc;
        for (int i = 0, ij = 0; i < fdim; ++i) {
            for (int j = 0; j < fdim; ++j, ++ij) {
                if (i == j) {
                    rho[ij] = (i == occrand) ? phys::math::iu : phys::math::iz;
                } else if (i == occrand || j == occrand) {
                    rho[ij] *= gamma_ou;
                } else {
                    rho[ij] *= gamma_uu;
                }
            }
        }
        return 0;
    }

    static int c_window(kids_complex* c, int iocc, int type, int fdim) {
        switch (type) {
            case SQCPolicy::TRI: {
                kids_real tmp2[2];
                Kernel_Random::rand_uniform(tmp2, 2);
                while (tmp2[0] + tmp2[1] > 1.0f) Kernel_Random::rand_uniform(tmp2, 2);
                c[iocc] = tmp2[0];
                tmp2[1] = 1.0f - std::real(c[iocc]);
                for (int i = 0; i < fdim; ++i) {
                    if (i != iocc) {
                        Kernel_Random::rand_uniform(tmp2, 1);
                        c[i] = tmp2[0] * tmp2[1];
                    }
                }
                c[iocc] += 1.0e0;
                break;
            }
            case SQCPolicy::SPX: {
                c_sphere(c, Dimension::F);
                for (int i = 0; i < Dimension::F; ++i) c[i] = std::abs(c[i] * c[i]);
                c[iocc] += 1.0e0;
                break;
            }
            case SQCPolicy::BIG: {
                kids_complex* cadd1 = new kids_complex[Dimension::Fadd1];
                c_sphere(cadd1, Dimension::Fadd1);
                for (int i = 0; i < Dimension::F; ++i) c[i] = std::abs(cadd1[i] * cadd1[i]);
                c[iocc] += 1.0e0;
                delete[] cadd1;
                break;
            }
            case SQCPolicy::SQR: {
                const kids_real gm0 = gamma_wigner(2.0f);
                for (int i = 0; i < fdim; ++i) {
                    kids_real randu;
                    Kernel_Random::rand_uniform(&randu);
                    c[i] = 2.0 * gm0 * randu;
                }
                c[iocc] += 1.0e0;
                break;
            }
        }
        for (int i = 0; i < fdim; ++i) {
            kids_real randu;
            Kernel_Random::rand_uniform(&randu);
            randu *= phys::math::twopi;
            c[i] = sqrt(c[i]);
            c[i] *= (cos(randu) + phys::math::im * sin(randu));
        }
        return 0;
    };

    static int ker_binning(kids_complex* ker, kids_complex* rho, int sqc_type) {
        const kids_real gm0 = gamma_wigner(2.0f), gm1 = 1 + gm0, gmh = 0.5f + gm0;

        // set all elements to 1
        for (int i = 0; i < Dimension::FF; ++i) ker[i] = rho[i] / std::abs(rho[i]);

        // for ker[ij], loop k-index to find if ker[ij] should be set to 0
        for (int i = 0, ij = 0; i < Dimension::F; ++i) {
            for (int j = 0; j < Dimension::F; ++j, ++ij) {
                for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1) {
                    double vk      = std::abs(rho[kk]);
                    bool   Outlier = false;
                    switch (sqc_type) {
                        case SQCPolicy::TRI:
                        case SQCPolicy::SPX:
                        case SQCPolicy::BIG:
                            Outlier = (i == j) ? ((k != i && vk > 1) || (k == i && vk < 1))
                                               : ((k != i && k != j && vk > 1) || ((k == i || k == j) && vk < 0.5f));

                            if (i != j) Outlier = false;
                            break;
                        case SQCPolicy::SQR:  // @bug?
                            Outlier = (i == j) ? ((k != i && std::abs(vk - gm0) < gm0) ||
                                                  (k == i && std::abs(vk - gm1) < gm0))
                                               : ((k != i && std::abs(vk - gm0) > gm0) ||  //
                                                  (k == i && std::abs(vk - gmh) > gm0) ||  //
                                                  (k == j && std::abs(vk - gmh) > gm0));
                            break;
                    }
                    if (Outlier) {
                        ker[ij] = phys::math::iz;
                        break;
                    }
                }
            }
        }
        return 0;
    };
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_Utils_H

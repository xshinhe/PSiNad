#include "Kernel_Elec_SQC.h"

#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Elec.h"
#include "../kernels/Kernel_Elec_CMM.h"
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

int Kernel_Elec_SQC::c_window(num_complex *c, int iocc, int type, int fdim) {
    switch (type) {
        case SQCPolicy::TRI: {
            num_real tmp2[2];
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
        case SQCPolicy::SQR: {
            const num_real gm0 = Kernel_Elec_CMM::gamma_wigner(2.0f);
            for (int i = 0; i < fdim; ++i) {
                num_real randu;
                Kernel_Random::rand_uniform(&randu);
                c[i] = 2.0 * gm0 * randu;
            }
            c[iocc] += 1.0e0;
            break;
        }
    }
    for (int i = 0; i < fdim; ++i) {
        num_real randu;
        Kernel_Random::rand_uniform(&randu);
        randu *= phys::math::twopi;
        c[i] *= (cos(randu) + phys::math::im * sin(randu));
    }
    return 0;
};

void Kernel_Elec_SQC::read_param_impl(Param *PM) {
    sqc_type = SQCPolicy::_from(PM->get<std::string>("sqc_flag", LOC(), "TRI"));
    switch (sqc_type) {
        case SQCPolicy::TRI:
            gamma = PM->get<double>("gamma", LOC(), 1.0 / 3);
            break;
        case SQCPolicy::SQR:
            gamma = PM->get<double>("gamma", LOC(), Kernel_Elec_CMM::gamma_wigner(2));
            break;
    }
    use_cv = PM->get<bool>("use_cv", LOC(), false);
}

void Kernel_Elec_SQC::init_calc_impl(int stat) {
    Kernel_Elec::w[0] = num_complex(1);
    c_window(Kernel_Elec::c, Kernel_Elec::occ0, sqc_type, Dimension::F);  // non-standard c

    *Kernel_Elec::occ_nuc = Kernel_Elec::occ0;                                          // useless
    Kernel_Elec::ker_from_c(Kernel_Elec::rho_ele, Kernel_Elec::c, 1, 0, Dimension::F);  // single-rank
    Kernel_Elec::ker_from_rho(Kernel_Elec::rho_nuc, Kernel_Elec::rho_ele, 1, gamma, Dimension::F);
    if (use_cv) {
        for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
            Kernel_Elec::rho_nuc[ii] = (i == Kernel_Elec::occ0) ? phys::math::iu : phys::math::iz;
        }
    }

    for (int i = 0, ik = 0; i < Dimension::F; ++i) {
        for (int k = 0; k < Dimension::F; ++k, ++ik) {
            Kernel_Elec::K1[ik] = (i == Kernel_Elec::occ0 && k == Kernel_Elec::occ0) ? phys::math::iu : phys::math::iz;
        }
    }

    exec_kernel(stat);
}

int Kernel_Elec_SQC::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::FF; ++i)
        Kernel_Elec::K2[i] = Kernel_Elec::rho_ele[i] / std::abs(Kernel_Elec::rho_ele[i]);

    // then set zeros (quantize to zero by window function)
    switch (sqc_type) {
        case SQCPolicy::TRI: {
            for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                for (int j = 0; j < Dimension::F; ++j, ++ij) {
                    for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1) {
                        double vk = std::abs(Kernel_Elec::rho_ele[kk]);
                        if ((i == j && ((k != i && vk > 1) || (k == i && vk < 1))) ||  // diagonal
                            (i != j && ((k != i && k != j && vk > 1) || ((k == i || k == j) && vk < 0.5f)))) {
                            Kernel_Elec::K2[ij] = phys::math::iz;
                            break;
                        }
                    }
                }
            }
            break;
        }
        case SQCPolicy::SQR: {
            const num_real gm0 = Kernel_Elec_CMM::gamma_wigner(2.0f), gm1 = 1 + gm0, gmh = 0.5f + gm0;
            for (int i = 0, ij = 0; i < Dimension::F; ++i) {
                for (int j = 0; j < Dimension::F; ++j, ++ij) {
                    for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1) {
                        double vk = std::abs(Kernel_Elec::rho_ele[kk]);
                        if ((i == j && ((k != i && std::abs(vk - gm0) > gm0) ||
                                        (k == i && std::abs(vk - gm1) > gm0))) ||  // diagonal
                            (i != j && ((k != i && std::abs(vk - gm0) > gm0) ||    //
                                        (k == i && std::abs(vk - gmh) > gm0) ||    //
                                        (k == j && std::abs(vk - gmh) > gm0)))) {
                            Kernel_Elec::K2[ij] = phys::math::iz;
                            break;
                        }
                    }
                }
            }
            break;
        }
    }
    return 0;
};



};  // namespace PROJECT_NS

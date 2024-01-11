#include "Kernel_Elec_SQC.h"

#include "../core/linalg.h"
#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Elec.h"
#include "../kernels/Kernel_Elec_CMM.h"
#include "../kernels/Kernel_Random.h"
#include "../kernels/Kernel_Representation.h"

namespace PROJECT_NS {

int Kernel_Elec_SQC::c_window(num_complex* c, int iocc, int type, int fdim) {
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
        case SQCPolicy::SPX: {
            Kernel_Elec_CMM::c_sphere(c, Dimension::F);
            for (int i = 0; i < Dimension::F; ++i) c[i] = abs(c[i] * c[i]);
            c[iocc] += 1.0e0;
            break;
        }
        case SQCPolicy::BIG: {
            num_complex* cadd1 = new num_complex[Dimension::Fadd1];
            Kernel_Elec_CMM::c_sphere(cadd1, Dimension::Fadd1);
            for (int i = 0; i < Dimension::F; ++i) c[i] = abs(cadd1[i] * cadd1[i]);
            c[iocc] += 1.0e0;
            delete[] cadd1;
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
        c[i] = sqrt(c[i]);
        c[i] *= (cos(randu) + phys::math::im * sin(randu));
    }
    return 0;
};

int Kernel_Elec_SQC::ker_binning(num_complex* ker, num_complex* rho, int sqc_type) {
    const num_real gm0 = Kernel_Elec_CMM::gamma_wigner(2.0f), gm1 = 1 + gm0, gmh = 0.5f + gm0;

    // set all elements to 1
    for (int i = 0; i < Dimension::FF; ++i) ker[i] = rho[i] / std::abs(rho[i]);

    // for ker[ij], loop k-index to find if ker[ij] should be set to 0
    for (int i = 0, ij = 0; i < Dimension::F; ++i) {
        for (int j = 0; j < Dimension::F; ++j, ++ij) {
            for (int k = 0, kk = 0; k < Dimension::F; ++k, kk += Dimension::Fadd1) {
                double vk    = std::abs(rho[kk]);
                bool Outlier = false;
                switch (sqc_type) {
                    case SQCPolicy::TRI:
                    case SQCPolicy::SPX:
                    case SQCPolicy::BIG:
                        Outlier = (i == j) ? ((k != i && vk > 1) || (k == i && vk < 1))
                                           : ((k != i && k != j && vk > 1) || ((k == i || k == j) && vk < 0.5f));

                        if (i != j) Outlier = false;

                        // if (abs(rho[i * Dimension::Fadd1] - rho[j * Dimension::Fadd1]) > 1.0f) Outlier = true;

                        // Outlier =
                        //     (i == j) ? ((k != i && vk > 1) || (k == i && vk < 1)) : (((k == i || k == j) && vk
                        //     > 1.0f));
                        break;
                    case SQCPolicy::SQR:  // @bug
                        Outlier = (i == j)
                                      ? ((k != i && std::abs(vk - gm0) < gm0) || (k == i && std::abs(vk - gm1) < gm0))
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

void Kernel_Elec_SQC::read_param_impl(Param* PM) {
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

void Kernel_Elec_SQC::init_data_impl(DataSet* DS) {
    sqcw  = _DataSet->reg<num_real>("integrator.sqcw", Dimension::FF);
    sqcw0 = _DataSet->reg<num_real>("integrator.sqcw0", Dimension::FF);
    sqcwh = _DataSet->reg<num_real>("integrator.sqcwh", Dimension::FF);
};

void Kernel_Elec_SQC::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* w       = Kernel_Elec::w + iP;
        num_complex* c       = Kernel_Elec::c + iP * Dimension::F;
        num_complex* rho_ele = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_nuc = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* U       = Kernel_Elec::U + iP * Dimension::FF;
        int* occ_nuc         = Kernel_Elec::occ_nuc + iP;

        /////////////////////////////////////////////////////////////////

        w[0]     = 1.0e0;                                         ///< initial measure
        *occ_nuc = Kernel_Elec::occ0;                             ///< initial occupation
        c_window(c, *occ_nuc, sqc_type, Dimension::F);            ///< initial c: non-standard c
        Kernel_Elec::ker_from_c(rho_ele, c, 1, 0, Dimension::F);  ///< initial rho_ele = cc^
        Kernel_Elec::ker_from_rho(rho_nuc, rho_ele, 1, gamma, Dimension::F, use_cv, *occ_nuc);  ///< initial rho_nuc
        ARRAY_EYE(U, Dimension::F);                                                             ///< initial propagator
    }
    Kernel_Elec::c_init       = _DataSet->set("init.c", Kernel_Elec::c, Dimension::PF);
    Kernel_Elec::rho_ele_init = _DataSet->set("init.rho_ele", Kernel_Elec::rho_ele, Dimension::PFF);
    Kernel_Elec::rho_nuc_init = _DataSet->set("init.rho_nuc", Kernel_Elec::rho_nuc, Dimension::PFF);
    Kernel_Elec::T_init       = _DataSet->set("init.T", Kernel_Elec::T, Dimension::PFF);
    exec_kernel(stat);
}

int Kernel_Elec_SQC::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_complex* U            = Kernel_Elec::U + iP * Dimension::FF;
        num_complex* rho_ele      = Kernel_Elec::rho_ele + iP * Dimension::FF;
        num_complex* rho_ele_init = Kernel_Elec::rho_ele_init + iP * Dimension::FF;
        num_complex* rho_nuc      = Kernel_Elec::rho_nuc + iP * Dimension::FF;
        num_complex* rho_nuc_init = Kernel_Elec::rho_nuc_init + iP * Dimension::FF;
        num_complex* K1           = Kernel_Elec::K1 + iP * Dimension::FF;
        num_complex* K2           = Kernel_Elec::K2 + iP * Dimension::FF;
        num_real* T               = Kernel_Elec::T + iP * Dimension::FF;
        num_real* T_init          = Kernel_Elec::T_init + iP * Dimension::FF;

        // 1) transform from inp_repr => ele_repr
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_ele[ik] = rho_ele_init[ik];
        for (int ik = 0; ik < Dimension::FF; ++ik) rho_nuc[ik] = rho_nuc_init[ik];
        Kernel_Representation::transform(rho_ele, T_init, Dimension::F,         //
                                         Kernel_Representation::inp_repr_type,  //
                                         Kernel_Representation::ele_repr_type,  //
                                         SpacePolicy::L);
        // Kernel_Representation::transform(rho_nuc, T_init, Dimension::F,         //
        //                                  Kernel_Representation::inp_repr_type,  //
        //                                  Kernel_Representation::ele_repr_type,  //
        //                                  SpacePolicy::L); // @debug

        // 2) propagte along ele_repr
        ARRAY_MATMUL3_TRANS2(rho_ele, U, rho_ele, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F);
        // ARRAY_MATMUL3_TRANS2(rho_nuc, U, rho_nuc, U, Dimension::F, Dimension::F, Dimension::F, Dimension::F); //
        // @debug

        // 3) transform back from ele_repr => inp_repr
        Kernel_Representation::transform(rho_ele, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);
        Kernel_Representation::transform(rho_nuc, T, Dimension::F,              //
                                         Kernel_Representation::ele_repr_type,  //
                                         Kernel_Representation::inp_repr_type,  //
                                         SpacePolicy::L);

        ker_binning(K1, rho_ele, sqc_type);
        for (int i = 0, ik = 0; i < Dimension::F; ++i) {
            for (int k = 0; k < Dimension::F; ++k, ++ik) {
                double radius = abs(2.0e0 - rho_ele[i * Dimension::Fadd1] - rho_ele[k * Dimension::Fadd1]);
                radius        = sqrt(radius * radius + 0.0025e0);
                switch (sqc_type) {
                    case SQCPolicy::SPX: {
                        sqcw[ik] = pow(radius, 3 - Dimension::F);
                        break;
                    }
                    default:
                    case SQCPolicy::BIG: {
                        sqcw[ik] = pow(radius, 2 - Dimension::F);
                        break;
                    }
                }
            }
        }
        for (int ik = 0; ik < Dimension::FF; ++ik) K2[ik] = K1[ik];
    }
    return 0;
};



};  // namespace PROJECT_NS

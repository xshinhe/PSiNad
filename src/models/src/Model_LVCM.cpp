#include "psnd/Model_LVCM.h"

#include "psnd/Kernel_Random.h"
#include "psnd/debug_utils.h"
#include "psnd/hash_fnv1a.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

namespace PROJECT_NS {

const std::string Model_LVCM::getName() { return "Model_LVCM"; }

int Model_LVCM::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_LVCM::setInputParam_impl(std::shared_ptr<Param> PM) {
    lvcm_type = LVCMPolicy::_from(_param->get_string({"model.lvcm_flag"}, LOC(), "PYR3"));
}

void Model_LVCM::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    Hsys = DS->def(DATA::model::Hsys);
    w    = DS->def(DATA::model::w);
    Kmat = DS->def(DATA::model::Kmat);
    Qmat = DS->def(DATA::model::Qmat);
    Tmod = DS->def(DATA::model::Tmod);
    memset(Hsys.data(), 0, Dimension::FF * sizeof(psnd_real));
    switch (lvcm_type) {
        case LVCMPolicy::PYR3: {
            psnd_assert(Dimension::N == 3, "Dimension Error");
            double H_unit = phys::au_2_ev;

            N_mode = 2;
            // parameter for PYR3
            double  E_data_PYR3[2]      = {3.94, 4.84};
            double  w_data_PYR3[3]      = {0.074, 0.126, 0.118};
            double  kcoeff_data_PYR3[4] = {-0.105, 0.149, 0.037, -0.254};
            double  lcoeff_data_PYR3[4] = {0.000, 0.262, 0.262, 0.000};
            double *E_data              = E_data_PYR3;
            double *w_data              = w_data_PYR3;
            double *kcoeff_data         = kcoeff_data_PYR3;
            double *lcoeff_data         = lcoeff_data_PYR3;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {  //
                Hsys[ii] = E_data[i] / H_unit;
            }
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, idxkcoeff = 0, idxlcoeff = 0; j < Dimension::N; ++j) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik) {
                        Qmat[j * Dimension::FF + ik] =                                      //
                            (j < N_mode) ? ((i == k) ? (kcoeff_data[idxkcoeff++]) : 0.0e0)  //
                                         : lcoeff_data[idxlcoeff++];
                        Qmat[j * Dimension::FF + ik] /= H_unit;
                        Qmat[j * Dimension::FF + ik] *= std::sqrt(w[j]);
                    }
                }
            }
            break;
        }
        case LVCMPolicy::PYR24: {
            psnd_assert(Dimension::N == 24, "Dimension Error");
            double H_unit = phys::au_2_ev;

            N_mode = 23;
            // parameter for PYR24
            double E_data_PYR24[2]       = {-0.4617, 0.4617};
            double w_data_PYR24[24]      = {0.0740, 0.1273, 0.1568, 0.1347, 0.3431, 0.1157, 0.3242, 0.3621,
                                            0.2673, 0.3052, 0.0968, 0.0589, 0.0400, 0.1726, 0.2863, 0.2484,
                                            0.1536, 0.2105, 0.0778, 0.2294, 0.1915, 0.4000, 0.3810, 0.0936};
            double kcoeff_data_PYR24[46] = {-0.0964, 0.1194,  0.0470, 0.2012,  0.1594, 0.0484,  0.0308, -0.0308,
                                            0.0782,  -0.0782, 0.0261, -0.0261, 0.0717, -0.0717, 0.0780, -0.0780,
                                            0.0560,  -0.0560, 0.0625, -0.0625, 0.0188, -0.0188, 0.0112, -0.0112,
                                            0.0069,  -0.0069, 0.0265, -0.0265, 0.0433, -0.0433, 0.0361, -0.0361,
                                            0.0210,  -0.0210, 0.0281, -0.0281, 0.0102, -0.0102, 0.0284, -0.0284,
                                            0.0196,  -0.0196, 0.0306, -0.0306, 0.0269, -0.0269};
            double lcoeff_data_PYR24[4]  = {0.000, 0.1825, 0.1825, 0.000};
            double *E_data               = E_data_PYR24;
            double *w_data               = w_data_PYR24;
            double *kcoeff_data          = kcoeff_data_PYR24;
            double *lcoeff_data          = lcoeff_data_PYR24;

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {  //
                Hsys[ii] = E_data[i] / H_unit;
            }
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, idxkcoeff = 0, idxlcoeff = 0; j < Dimension::N; ++j) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik) {
                        Qmat[j * Dimension::FF + ik] =                                      //
                            (j < N_mode) ? ((i == k) ? (kcoeff_data[idxkcoeff++]) : 0.0e0)  //
                                         : lcoeff_data[idxlcoeff++];
                        Qmat[j * Dimension::FF + ik] /= H_unit;
                        Qmat[j * Dimension::FF + ik] *= std::sqrt(w[j]);
                    }
                }
            }
            break;
        }
        case LVCMPolicy::BUTA5: {
            psnd_assert(Dimension::N == 5, "Dimension Error");
            double H_unit = phys::au_2_ev;

            N_mode = 4;
            // parameter for BUTA5
            double  E_data_BUTA5[2]      = {9.41165, 9.95575};
            double  w_data_BUTA5[5]      = {0.1089, 0.1773, 0.2578, 0.3713, 0.0912};
            double  kcoeff_data_BUTA5[8] = {-0.0531, -0.0594, 0.0115,  0.0100,  //
                                            -0.1628, 0.3422,  -0.0403, 0.0321};
            double  lcoeff_data_BUTA5[4] = {0.000, 0.2880, 0.2880, 0.000};
            double *E_data               = E_data_BUTA5;
            double *w_data               = w_data_BUTA5;
            double *kcoeff_data          = kcoeff_data_BUTA5;
            double *lcoeff_data          = lcoeff_data_BUTA5;

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {  //
                Hsys[ii] = E_data[i] / H_unit;
            }
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, idxkcoeff = 0, idxlcoeff = 0; j < Dimension::N; ++j) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik) {
                        Qmat[j * Dimension::FF + ik] =                                      //
                            (j < N_mode) ? ((i == k) ? (kcoeff_data[idxkcoeff++]) : 0.0e0)  //
                                         : lcoeff_data[idxlcoeff++];
                        Qmat[j * Dimension::FF + ik] /= H_unit;
                        Qmat[j * Dimension::FF + ik] *= std::sqrt(w[j]);
                    }
                }
            }
            break;
        }
        case LVCMPolicy::CRC2:
        case LVCMPolicy::CRC5: {
            psnd_assert(Dimension::N <= 5, "Dimension Error");
            double H_unit = phys::au_2_ev;

            N_mode = 0;

            const double kappa1   = -0.0169;
            const double kappa3   = -0.0272;
            const double lambda1a = 0.0328;
            const double lambda1b = -0.0978;
            const double lambda2a = 0.0095;
            const double lambda2b = 0.1014;

            double E_data_CRC5[3]       = {0.0424, 0.0424, 0.4344};
            double w_data_CRC5[5]       = {0.0129, 0.0129, 0.0342, 0.0561, 0.0561};
            double kcoeff_data_CRC5[1]  = {0.0};
            double lcoeff_data_CRC5[45] = {
                // data
                0.0,      lambda1a, 0.0,      // x0
                lambda1a,  0.0,     lambda1b,  //
                0.0,      lambda1b, 0.0,      //
                -lambda1a, 0.0,     lambda1b,  // x1
                0.0,      lambda1a, 0.0,      //
                lambda1b,  0.0,     0.0,      //
                kappa1,    0.0,     0.0,      // x2
                0.0,      kappa1,   0.0,      //
                0.0,      0.0,     kappa3,    //
                -lambda2a, 0.0,     lambda2b,  // x3
                0.0,      lambda2a, 0.0,      //
                lambda2b,  0.0,     0.0,      //
                0.0,      lambda2a, 0.0,      // x4
                lambda2a,  0.0,     lambda2b,  //
                0.0,      lambda2b, 0.0       //
            };

            double *E_data      = E_data_CRC5;
            double *w_data      = w_data_CRC5;
            double *kcoeff_data = kcoeff_data_CRC5;
            double *lcoeff_data = lcoeff_data_CRC5;

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {  //
                Hsys[ii] = E_data[i] / H_unit;
            }
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, idxkcoeff = 0, idxlcoeff = 0; j < Dimension::N; ++j) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik) {
                        Qmat[j * Dimension::FF + ik] =                                      //
                            (j < N_mode) ? ((i == k) ? (kcoeff_data[idxkcoeff++]) : 0.0e0)  //
                                         : lcoeff_data[idxlcoeff++];
                        Qmat[j * Dimension::FF + ik] /= H_unit;
                        Qmat[j * Dimension::FF + ik] *= std::sqrt(w[j]);
                    }
                }
            }
            break;
        }
        case LVCMPolicy::CED2:
        case LVCMPolicy::CED3: {
            N_mode                  = 0;
            const double lightspeed = 137.03599907444;
            const double epsilon0   = 0.25 / phys::math::pi;

            double E_data_CED2[2]  = {-0.6738, -0.2798};
            double mu_data_CED2[4] = {0.000, +1.034,  //
                                      +1.034, 0.000};

            double E_data_CED3[3]  = {-0.6738, -0.2798, -0.1547};
            double mu_data_CED3[9] = {0.000,  +1.034, 0.000,   //
                                      +1.034, 0.000,  -2.536,  //
                                      0.000,  -2.536, 0.000};

            double Lcav = 2.362e5;
            double Rcav = Lcav / 2;

            double *E_data;
            double *mu_data;
            switch (lvcm_type) {
                case LVCMPolicy::CED2:
                    psnd_assert(Dimension::F == 2, "Dimension Error");

                    E_data  = E_data_CED2;
                    mu_data = mu_data_CED2;
                    break;
                case LVCMPolicy::CED3:
                    psnd_assert(Dimension::F == 3, "Dimension Error");

                    E_data  = E_data_CED3;
                    mu_data = mu_data_CED3;
                    break;
            }

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) { Hsys[ii] = E_data[i]; }
            for (int j = 0, jik = 0; j < Dimension::N; ++j) {
                w[j] = (2 * j + 1) * lightspeed * phys::math::pi / Lcav;
                for (int ik = 0; ik < Dimension::FF; ++ik, ++jik) {
                    Qmat[jik] = sqrt(2.0 / (epsilon0 * Lcav)) * sin((2 * j + 1) * phys::math::pi * Rcav / Lcav) *
                                w[j] * mu_data[ik];
                }
            }
            break;
        }
        case LVCMPolicy::PYR2CED: {
            N_mode        = 1;
            double H_unit = phys::au_2_ev;

            double gcoup = _param->get_real({"model.gcoup"}, LOC(), phys::energy_d, 0.24f / H_unit);
            double wcav  = _param->get_real({"model.wcav"}, LOC(), phys::energy_d, 0.62f / H_unit);

            // parameter for PYR2
            double  E_data_PYR2[2]      = {3.94, 4.84};
            double  w_data_PYR2[3]      = {0.074, 0.118};
            double  kcoeff_data_PYR2[4] = {-0.105, 0.149};
            double  lcoeff_data_PYR2[4] = {0.000, 0.262, 0.262, 0.000};
            double *E_data              = E_data_PYR2;
            double *w_data              = w_data_PYR2;
            double *kcoeff_data         = kcoeff_data_PYR2;
            double *lcoeff_data         = lcoeff_data_PYR2;

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                Hsys[ii] = i / (Dimension::F / 2) * wcav + E_data[i % 2] / H_unit;
            }
            Hsys[0 * 4 + 3] = gcoup;
            Hsys[1 * 4 + 2] = gcoup;
            Hsys[2 * 4 + 1] = gcoup;
            Hsys[3 * 4 + 0] = gcoup;
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, ji = 0; j < N_mode; ++j) {
                for (int i = 0; i < Dimension::F; ++i, ++ji) {
                    kcoeff[ji] = (kcoeff_data[j * 2 + i % 2] / H_unit) * sqrt(w[j]);
                }
            }
            for (int j = 0, jik = 0; j < Dimension::N; ++j) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik, ++jik) {
                        if (j < N_mode) {
                            Qmat[jik] = (i == k) ? (kcoeff_data[j * 2 + i % 2] / H_unit) * sqrt(w[j]) : 0.0e0;
                        } else {
                            Qmat[jik] =
                                (i / (Dimension::F / 2) == k / (Dimension::F / 2))  //
                                    ? (lcoeff_data[(j - N_mode) * 4 + (i % 2) * 2 + k % 2] / H_unit) * sqrt(w[j])
                                    : 0.0e0;
                        }
                    }
                }
            }
            break;
        }
        case LVCMPolicy::BEN5: {
            break;
        }
        case LVCMPolicy::Read: {
            std::string   lvcm_file = _param->get_string({"model.lvcm_file"}, LOC(), "lvcm.dat");
            std::ifstream ifs(lvcm_file);
            psnd_assert(ifs.is_open(), "File not found: " + lvcm_file);   // check if file exists
            std::string   H_unit_str;
            std::string   firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit
            double H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));

            // read E
            int         dsize;
            std::string flag;
            double      val;
            ifs >> flag >> dsize;
            psnd_assert(dsize == Dimension::F, "Dimension Error");
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1)
                if (ifs >> val) Hsys[ii] = val / H_unit;

            // read w
            ifs >> flag >> dsize;
            psnd_assert(dsize == Dimension::N, "Dimension Error");
            for (int i = 0, ii = 0; i < Dimension::N; ++i)
                if (ifs >> val) w[i] = val / H_unit;

            N_mode = 0;
            // read kcoeff & lcoeff
            for (int j = 0; j < Dimension::N; ++j) {
                ifs >> flag;
                if (flag == "K") {
                    for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                        if (ifs >> val) {
                            Qmat[j * Dimension::FF + ii] = val / H_unit;
                            Qmat[j * Dimension::FF + ii] *= std::sqrt(w[j]);
                        }
                    }
                    N_mode++;
                } else if (flag == "L") {
                    for (int ik = 0; ik < Dimension::FF; ++ik) {
                        if (ifs >> val) {
                            Qmat[j * Dimension::FF + ik] = val / H_unit;
                            Qmat[j * Dimension::FF + ik] *= std::sqrt(w[j]);
                        }
                    }
                } else {  // not modified with sqrt(w)
                    for (int ik = 0; ik < Dimension::FF; ++ik) {
                        if (ifs >> val) { Qmat[j * Dimension::FF + ik] = val / H_unit; }
                    }
                }
            }
            ifs.close();
        }
    }

    /// 2) init Bath sub-kernel (declaration & call)
    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);
    switch (lvcm_type) {
        case LVCMPolicy::CRC2:
        case LVCMPolicy::CRC5: {
            double reqb[5] = {0.0, 14.3514, -9.9699, -7.0189, 0.0};
            double alpw[5] = {0.4501, 0.4286, 0.6204, 0.4535, 0.5539};
            for (int j = 0; j < Dimension::N; ++j) {
                x0[j]      = reqb[j] / sqrt(w[j]);
                p0[j]      = 0.0;
                x_sigma[j] = alpw[j] / sqrt(w[j]);
                p_sigma[j] = 0.5 * sqrt(w[j]) / alpw[j];
            }
            break;
        }
        default: {
            for (int j = 0; j < Dimension::N; ++j) {
                x0[j]      = 0.0;
                p0[j]      = 0.0;
                x_sigma[j] = sqrt(0.5 / w[j]);
                p_sigma[j] = sqrt(0.5 * w[j]);
            }
            break;
        }
    }

    for (int j = 0, jk = 0; j < Dimension::N; ++j) {
        for (int k = 0; k < Dimension::N; ++k, ++jk) {
            Kmat[jk] = (j == k) ? w[j] * w[j] : 0.0e0;
            Tmod[jk] = (j == k) ? 1.0e0 : 0.0e0;
        }
    }

    // model field
    mass = DS->def(DATA::model::mass);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0;
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

Status &Model_LVCM::initializeKernel_impl(Status &stat) { return stat; }

Status &Model_LVCM::executeKernel_impl(Status &stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto x    = this->x.subspan(iP * Dimension::N, Dimension::N);
        auto vpes = this->vpes.subspan(iP, 1);
        auto grad = this->grad.subspan(iP * Dimension::N, Dimension::N);
        auto hess = this->hess.subspan(iP * Dimension::NN, Dimension::NN);
        auto V    = this->V.subspan(iP * Dimension::FF, Dimension::FF);
        auto dV   = this->dV.subspan(iP * Dimension::NFF, Dimension::NFF);

        // calculate nuclear vpes and grad
        double term = 0.0;
        for (int j = 0; j < Dimension::N; ++j) {
            term += w[j] * w[j] * x[j] * x[j];
            grad[j] = w[j] * w[j] * x[j];
        }
        vpes[0] = 0.5 * term;

        // electronic pes
        memset(V.data(), 0, Dimension::FF * sizeof(psnd_real));
        for (int ik = 0; ik < Dimension::FF; ++ik) V[ik] = Hsys[ik];
        // ARRAY_SHOW(V, Dimension::F, Dimension::F);

        for (int j = 0, jFF = 0; j < N_mode; ++j, jFF += Dimension::FF) {
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) {
                V[ii] += Qmat[jFF + ii] * x[j];  // fast
            };
        }
        for (int j = N_mode, jik = N_mode * Dimension::FF; j < Dimension::N; ++j) {
            for (int ik = 0; ik < Dimension::FF; ++ik, ++jik) {
                V[ik] += Qmat[jik] * x[j];  // slow
            };
        }

        // ARRAY_SHOW(V, Dimension::F, Dimension::F);

        if (count_exec == 0) {
            for (int i = 0; i < Dimension::NFF; ++i) dV[i] = Qmat[i];
            // ddV = 0;
        }
    }
    return stat;
}


};  // namespace PROJECT_NS

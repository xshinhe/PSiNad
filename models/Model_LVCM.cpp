#include "Model_LVCM.h"

#include "../kernels/Kernel_Declare.h"
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

void Model_LVCM::read_param_impl(Param *PM) {
    lvcm_type = LVCMPolicy::_from(_Param->get<std::string>("lvcm_flag", LOC(), "PYR3"));
}

void Model_LVCM::init_data_impl(DataSet *DS) {
    Hsys = DS->reg<num_real>("model.Hsys", Dimension::FF);
    w    = DS->reg<num_real>("model.w", Dimension::N);
    memset(Hsys, 0, Dimension::FF * sizeof(num_real));
    switch (lvcm_type) {
        case LVCMPolicy::PYR3:
        case LVCMPolicy::PYR24: {
            double H_unit = phys::au_2_ev;

            // parameter for PYR3
            double E_data_PYR3[2]      = {3.94f, 4.84f};
            double w_data_PYR3[3]      = {0.074f, 0.126f, 0.118f};
            double kcoeff_data_PYR3[4] = {-0.105f, 0.149f, 0.037f, -0.254f};
            double lcoeff_data_PYR3[4] = {0.000f, 0.262f, 0.262f, 0.000f};

            // parameter for PYR24
            double E_data_PYR24[2]       = {-0.4617f, 0.4617f};
            double w_data_PYR24[24]      = {0.0740f, 0.1273f, 0.1568f, 0.1347f, 0.3431f, 0.1157f, 0.3242f, 0.3621f,
                                       0.2673f, 0.3052f, 0.0968f, 0.0589f, 0.0400f, 0.1726f, 0.2863f, 0.2484f,
                                       0.1536f, 0.2105f, 0.0778f, 0.2294f, 0.1915f, 0.4000f, 0.3810f, 0.0936f};
            double kcoeff_data_PYR24[46] = {-0.0964f, 0.1194f,  0.0470f, 0.2012f,  0.1594f, 0.0484f,  0.0308f, -0.0308f,
                                            0.0782f,  -0.0782f, 0.0261f, -0.0261f, 0.0717f, -0.0717f, 0.0780f, -0.0780f,
                                            0.0560f,  -0.0560f, 0.0625f, -0.0625f, 0.0188f, -0.0188f, 0.0112f, -0.0112f,
                                            0.0069f,  -0.0069f, 0.0265f, -0.0265f, 0.0433f, -0.0433f, 0.0361f, -0.0361f,
                                            0.0210f,  -0.0210f, 0.0281f, -0.0281f, 0.0102f, -0.0102f, 0.0284f, -0.0284f,
                                            0.0196f,  -0.0196f, 0.0306f, -0.0306f, 0.0269f, -0.0269f};
            double lcoeff_data_PYR24[4]  = {0.000f, 0.1825f, 0.1825f, 0.000f};


            double *E_data;
            double *w_data;
            double *kcoeff_data;
            double *lcoeff_data;
            switch (lvcm_type) {
                case LVCMPolicy::PYR3:
                    N_mode      = 2;
                    N_coup      = 1;
                    E_data      = E_data_PYR3;
                    w_data      = w_data_PYR3;
                    kcoeff_data = kcoeff_data_PYR3;
                    lcoeff_data = lcoeff_data_PYR3;
                    break;
                case LVCMPolicy::PYR24:
                    N_mode      = 23;
                    N_coup      = 1;
                    E_data      = E_data_PYR24;
                    w_data      = w_data_PYR24;
                    kcoeff_data = kcoeff_data_PYR24;
                    lcoeff_data = lcoeff_data_PYR24;
                    break;
            }

            kcoeff = DS->reg<num_real>("model.kcoeff", N_mode * Dimension::F);
            lcoeff = DS->reg<num_real>("model.lcoeff", N_coup * Dimension::FF);

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) { Hsys[ii] = E_data[i] / H_unit; }
            for (int j = 0; j < Dimension::N; ++j) w[j] = w_data[j] / H_unit;
            for (int j = 0, ji = 0; j < N_mode; ++j) {
                for (int i = 0; i < Dimension::F; ++i, ++ji) { kcoeff[ji] = (kcoeff_data[ji] / H_unit) * sqrt(w[j]); }
            }
            for (int j = N_mode, j0ik = 0; j < Dimension::N; ++j) {
                for (int ik = 0; ik < Dimension::FF; ++ik, ++j0ik) {
                    lcoeff[j0ik] = (lcoeff_data[j0ik] / H_unit) * sqrt(w[j]);
                }
            }
            break;
        }
        case LVCMPolicy::CED2:
        case LVCMPolicy::CED3: {
            const double lightspeed = 137.03599907444f;
            const double epsilon0   = 0.25f / phys::math::pi;

            double E_data_CED2[2]  = {-0.6738f, -0.2798f};
            double mu_data_CED2[4] = {0.000f, +1.034f,  //
                                      +1.034f, 0.000f};

            double E_data_CED3[3]  = {-0.6738f, -0.2798f, -0.1547f};
            double mu_data_CED3[9] = {0.000f,  +1.034f, 0.000f,   //
                                      +1.034f, 0.000f,  -2.536f,  //
                                      0.000f,  -2.536f, 0.000f};

            double Lcav = 2.362e5;
            double Rcav = Lcav / 2;

            double *E_data;
            double *mu_data;
            switch (lvcm_type) {
                case LVCMPolicy::CED2:
                    E_data  = E_data_CED2;
                    mu_data = mu_data_CED2;
                    break;
                case LVCMPolicy::CED3:
                    E_data  = E_data_CED3;
                    mu_data = mu_data_CED3;
                    break;
            }

            N_mode = 0;
            N_coup = Dimension::N;

            lcoeff = DS->reg<num_real>("model.lcoeff", N_coup * Dimension::FF);

            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1) { Hsys[ii] = E_data[i]; }
            for (int j = 0, jik = 0; j < Dimension::N; ++j) {
                w[j] = (2 * j + 1) * lightspeed * phys::math::pi / Lcav;
                for (int ik = 0; ik < Dimension::FF; ++ik, ++jik) {
                    lcoeff[jik] = sqrt(2.0f / (epsilon0 * Lcav)) * sin((2 * j + 1) * phys::math::pi * Rcav / Lcav) *
                                  w[j] * mu_data[ik];
                }
            }
            break;
        }
        case LVCMPolicy::PYR2CED: {
            N_mode = 1;
            N_coup = 1;

            double H_unit = phys::au_2_ev;

            double gcoup = _Param->get<double>("gcoup", LOC(), phys::energy_d, 0.24f / H_unit);
            double wcav  = _Param->get<double>("wcav", LOC(), phys::energy_d, 0.62f / H_unit);

            // parameter for PYR2
            double E_data_PYR2[2]      = {3.94f, 4.84f};
            double w_data_PYR2[3]      = {0.074f, 0.118f};
            double kcoeff_data_PYR2[4] = {-0.105f, 0.149f};
            double lcoeff_data_PYR2[4] = {0.000f, 0.262f, 0.262f, 0.000f};

            double *E_data      = E_data_PYR2;
            double *w_data      = w_data_PYR2;
            double *kcoeff_data = kcoeff_data_PYR2;
            double *lcoeff_data = lcoeff_data_PYR2;

            kcoeff = DS->reg<num_real>("model.kcoeff", N_mode * Dimension::F);
            lcoeff = DS->reg<num_real>("model.lcoeff", N_coup * Dimension::FF);

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
            for (int j = N_mode, j0 = 0, j0ik = 0; j < Dimension::N; ++j, ++j0) {
                for (int i = 0, ik = 0; i < Dimension::F; ++i) {
                    for (int k = 0; k < Dimension::F; ++k, ++ik, ++j0ik) {
                        if (i / (Dimension::F / 2) == k / (Dimension::F / 2))
                            lcoeff[j0ik] = (lcoeff_data[j0 * 4 + (i % 2) * 2 + k % 2] / H_unit) * sqrt(w[j]);
                    }
                }
            }
            // ARRAY_SHOW(Hsys, Dimension::F, Dimension::F);
            // ARRAY_SHOW(kcoeff, N_mode, Dimension::F);
            // ARRAY_SHOW(lcoeff, N_coup, Dimension::FF);
            break;
        }
        case LVCMPolicy::BEN5: {
            break;
        }
        case LVCMPolicy::Read: {
            std::string lvcm_readfile = _Param->get<std::string>("lvcm_readfile", LOC(), "lvcm.dat");
            std::ifstream ifs(lvcm_readfile);
            std::string H_unit_str;
            std::string firstline;
            getline(ifs, firstline);
            std::stringstream sstr(firstline);
            sstr >> H_unit_str;  ///< the firstline stores H's unit
            double H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));

            // read E
            int dsize;
            std::string flag;
            double val;
            ifs >> flag >> dsize;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1)
                if (ifs >> val) Hsys[ii] = val / H_unit;

            // read kcoeff & lcoeff
            for (int j = 0, ji = 0, j0ik = 0; j < Dimension::N; ++j) {
                ifs >> flag >> dsize;
                if (dsize == Dimension::F) {
                    if (ifs >> val) w[j] = val / H_unit;  // freq
                    for (int i = 0; i < Dimension::F; ++i, ++ji)
                        if (ifs >> val) kcoeff[ji] = val / H_unit;
                } else if (dsize = Dimension::FF) {
                    if (ifs >> val) w[j] = val / H_unit;  // freq
                    for (int ik = 0; ik < Dimension::FF; ++ik, ++j0ik)
                        if (ifs >> val) lcoeff[j0ik] = val / H_unit;
                } else {
                    exit(-1);
                }
            }
            ifs >> flag >> dsize;
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ii += Dimension::Fadd1)
                if (ifs >> val) Hsys[ii] = val / H_unit;

            for (int i = 0; i < Dimension::FF; ++i)
                if (ifs >> val) Hsys[i] = val / H_unit;
            ifs.close();
        }
    }

    /// 2) init Bath sub-kernel (declaration & call)
    x_sigma = DS->reg<double>("model.x_sigma", Dimension::N);
    p_sigma = DS->reg<double>("model.p_sigma", Dimension::N);
    for (int j = 0; j < Dimension::N; ++j) {
        x_sigma[j] = sqrt(0.5f / w[j]);
        p_sigma[j] = sqrt(0.5f * w[j]);
    }

    // model field
    mass = DS->reg<double>("model.mass", Dimension::N);
    for (int j = 0; j < Dimension::N; ++j) mass[j] = 1.0f;
    vpes = DS->reg<double>("model.vpes", Dimension::P);
    grad = DS->reg<double>("model.grad", Dimension::PN);
    hess = DS->reg<double>("model.hess", Dimension::PNN);
    V    = DS->reg<double>("model.V", Dimension::PFF);
    dV   = DS->reg<double>("model.dV", Dimension::PNFF);
    // ddV  = DS->reg<double>("model.ddV", Dimension::NNFF);

    // init & integrator
    x = DS->reg<double>("integrator.x", Dimension::PN);
    p = DS->reg<double>("integrator.p", Dimension::PN);
}

void Model_LVCM::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_real *x = this->x + iP * Dimension::N;
        num_real *p = this->p + iP * Dimension::N;

        Kernel_Random::rand_gaussian(x, Dimension::N);
        Kernel_Random::rand_gaussian(p, Dimension::N);
        for (int j = 0; j < Dimension::N; ++j) {
            x[j] = x[j] * x_sigma[j];
            p[j] = p[j] * p_sigma[j];
        }
    }

    _DataSet->set("init.x", x, Dimension::PN);
    _DataSet->set("init.p", p, Dimension::PN);
    exec_kernel(stat);
}

int Model_LVCM::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_real *x    = this->x + iP * Dimension::N;
        num_real *vpes = this->vpes + iP;
        num_real *grad = this->grad + iP * Dimension::N;
        num_real *hess = this->hess + iP * Dimension::NN;
        num_real *V    = this->V + iP * Dimension::FF;
        num_real *dV   = this->dV + iP * Dimension::NFF;
        // ARRAY_SHOW(x, 1, Dimension::N);

        // calculate nuclear vpes and grad
        double term = 0.0f;
        for (int j = 0; j < Dimension::N; ++j) {
            term += w[j] * w[j] * x[j] * x[j];
            grad[j] = w[j] * w[j] * x[j];
        }
        vpes[0] = 0.5 * term;

        // electronic pes
        memset(V, 0, Dimension::FF * sizeof(num_real));
        for (int ik = 0; ik < Dimension::FF; ++ik) V[ik] = Hsys[ik];
        // ARRAY_SHOW(V, Dimension::F, Dimension::F);

        for (int j = 0, ji = 0; j < N_mode; ++j) {
            for (int i = 0, ii = 0; i < Dimension::F; ++i, ++ji, ii += Dimension::Fadd1) {
                V[ii] += kcoeff[ji] * x[j];  //
            };
        }
        for (int j = N_mode, j0ik = 0; j < Dimension::N; ++j) {
            for (int ik = 0; ik < Dimension::FF; ++ik, ++j0ik) {
                V[ik] += lcoeff[j0ik] * x[j];  // slow
            };
        }

        // ARRAY_SHOW(V, Dimension::F, Dimension::F);

        if (count_exec == 0) {
            // N_mode
            for (int j = 0, ji = 0, jFF = 0; j < N_mode; ++j, jFF += Dimension::FF) {
                for (int i = 0, jii = jFF; i < Dimension::F; ++i, ++ji, jii += Dimension::Fadd1) {
                    dV[jii] = kcoeff[ji];  //
                }
            }
            // N_coup
            for (int j = N_mode, jFF = N_mode * Dimension::FF, j0FF = 0; j < Dimension::N;
                 ++j, jFF += Dimension::FF, j0FF += Dimension::FF) {
                for (int ik = 0, jik = jFF, j0ik = j0FF; ik < Dimension::FF; ++ik, ++jik, ++j0ik) {
                    dV[jik] = lcoeff[j0ik];
                }
            }
        }

        // ARRAY_SHOW(dV, Dimension::N, Dimension::FF);
        // exit(-1);
    }
    return 0;
}


};  // namespace PROJECT_NS

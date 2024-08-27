#include "kids/Model_HarmonicBath.h"

#include <algorithm>

#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_HarmonicBath::getName() { return "Model_HarmonicBath"; }

int Model_HarmonicBath::getType() const { return utils::hash(FUNCTION_NAME); }

double Model_HarmonicBath::J_Debye(double w) { return 2 * lambda * omegac * w / (w * w + omegac * omegac); }

double Model_HarmonicBath::J_Ohmic(double w) { return phys::math::pi * lambda / omegac * w * std::exp(-w / omegac); }

double Model_HarmonicBath::J(double w, double* w_arr, double* c_arr, int Nb) {
    if (w < 0) return -J(-w, w_arr, c_arr, Nb);
    if (w_arr == nullptr || c_arr == nullptr || Nb <= 0 || w / w_arr[Nb - 1] > 0.1)
        throw std::runtime_error("J(w) fitting error!");

    double a  = 10e0 * Nb * Nb / (w_arr[Nb - 1] * w_arr[Nb - 1]);
    double Jw = 0.0;
    for (int j = 0; j < Nb; ++j) {
        Jw += c_arr[j] * c_arr[j] / w_arr[j] * std::exp(-a * (w - w_arr[j]) * (w - w_arr[j]));
    }
    Jw *= 0.5e0 * phys::math::halfpi * std::sqrt(a / phys::math::pi);
    return Jw;
}

int Model_HarmonicBath::fun_Cw(kids_complex* Cw, double* w, int Nw, double* w_arr, double* c_arr, double beta, int Nb) {
    double dw    = std::min({1.0e-5, w_arr[0] / 10, w_arr[Nb - 1] / 1000});
    double C0_Re = (J(dw, w_arr, c_arr, Nb) - J(-dw, w_arr, c_arr, Nb)) / (2 * beta * dw);
    for (int i = 0; i < Nw; ++i) {
        // the real part
        Cw[i] = (w[i] == 0) ? C0_Re : J(w[i], w_arr, c_arr, Nb) * (1 + 1 / (exp(beta * w[i]) - 1.0f));
        // the imaginary part
        if (w[i] != 0) {
            double sum = 0.0e0, tmp = 1.0e0, Cp, Cm;
            double wp = w[i], wm = w[i];
            int    k = 0;
            while (std::abs(tmp) > 1.0e-6 || std::abs(tmp / sum) > 1.0e-8) {
                wp += dw, wm -= dw, k++;
                Cp  = (wp == 0) ? C0_Re : J(wp, w_arr, c_arr, Nb) * (1 + 1 / (exp(beta * wp) - 1.0f));
                Cm  = (wm == 0) ? C0_Re : J(wm, w_arr, c_arr, Nb) * (1 + 1 / (exp(beta * wm) - 1.0f));
                tmp = (Cp - Cm) / (double) k;
                sum += tmp;
            }
            Cw[i] += phys::math::im / phys::math::pi * (-sum);
        }
    }
    return 0;
}

void Model_HarmonicBath::setInputParam_impl(std::shared_ptr<Param> PM) {
    // size information
    int N = _param->get_int({"model.N"}, LOC());
    Nb    = _param->get_int({"model.Nb", "model.bath.Nb"}, LOC());
    nbath = _param->get_int({"model.nbath", "model.bath.nbath"}, LOC());
    kids_assert(Nb == Dimension::Nb, "Dimension Error");
    kids_assert(N == Dimension::N, "Dimension Error");

    is_classical    = _param->get_bool({"model.bath_classical", "model.bath.classical"}, LOC(), false);
    is_correlated   = _param->get_bool({"model.bath_correlated", "model.bath.correlated"}, LOC(), false);
    is_et_transform = _param->get_bool({"model.bath_et_transform", "model.bath.et_transform"}, LOC(), false);
    bath_type       = HarmonicBathPolicy::_from(_param->get_string({"model.bath_flag", "model.bath.flag"},  //
                                                             LOC(), "Debye"));
    strength_type   = StrengthPolicy::_from(_param->get_string({"model.strength_flag"}, LOC(), "Lambda"));
    omegac          = _param->get_real({"model.bath_omegac", "model.bath.omegac"}, LOC(), phys::energy_d, 1.0f);

    switch (strength_type) {
        case StrengthPolicy::Lambda: {
            double strength =
                _param->get_real({"model.bath_strength", "model.bath.strength"}, LOC(), phys::energy_d, 1.0f);
            lambda = strength;
            break;
        }
        case StrengthPolicy::Alpha: {
            double strength = _param->get_real({"model.bath_strength", "model.bath.strength"}, LOC(), 1.0f);
            lambda          = 0.5e0 * omegac * strength;
            break;
        }
        case StrengthPolicy::Eta: {
            double strength =
                _param->get_real({"model.bath_strength", "model.bath.strength"}, LOC(), phys::energy_d, 1.0f);
            lambda = 0.5e0 * strength;
            break;
        }
        case StrengthPolicy::Erg: {
            double strength =
                _param->get_real({"model.bath_strength", "model.bath.strength"}, LOC(), phys::energy_d, 1.0f);
            lambda = 0.25f * strength;
            break;
        }
    }
    double temperature =
        _param->get_real({"model.bath_temperature", "model.bath.temperature"}, LOC(), phys::temperature_d, 1.0f);
    beta = 1.0f / (phys::au::k * temperature);  // don't ignore k_Boltzman
}

void Model_HarmonicBath::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    w      = DS->def(DATA::model::w);
    Kmat   = DS->def(DATA::model::Kmat);
    Tmod   = DS->def(DATA::model::Tmod);
    omegas = DS->def(DATA::model::bath::omegas);
    coeffs = DS->def(DATA::model::bath::coeffs);
    if (nbath <= 0) {  // read from file
        std::string bath_file = _param->get_string({"model.bath_file", "model.bath.file"}, LOC(), "bath.dat");
        if (is_correlated) {
            try {
                std::ifstream ifs(bath_file);
                double        tmp;
                for (int j = 0; j < Dimension::NN; ++j) {
                    if (ifs >> tmp) Kmat[j] = tmp;
                }
                ifs.close();
            } catch (std::runtime_error& e) { throw kids_error("read Kmat.dat from bath_file fails"); }
            EigenSolve(w.data(), Tmod.data(), Kmat.data(), Dimension::N);
            for (int i = 0; i < Dimension::N; ++i) w[i] = std::sqrt(w[i]);
        } else {
            try {
                std::ifstream ifs(bath_file);
                double        tmp;
                for (int j = 0; j < Dimension::N; ++j) {
                    if (ifs >> tmp) w[j] = tmp;
                }
                ifs.close();
            } catch (std::runtime_error& e) { throw kids_error("read w.dat from bath_file fails"); }
            for (int j = 0, jk = 0; j < Dimension::N; ++j) {
                for (int k = 0; k < Dimension::N; ++k, ++jk) {
                    Kmat[jk] = (j == k) ? w[j] * w[j] : 0.0e0;
                    Tmod[jk] = (j == k) ? 1.0e0 : 0.0e0;
                }
            }
        }
    } else {  // built from spectrum function
        kids_assert(nbath >= 1, "spectrum function need nbath");
        kids_assert(!is_correlated, "spectrum function need diagonal un-correlated");
        switch (bath_type) {
            case HarmonicBathPolicy::Debye: {
                for (int j = 0; j < Dimension::Nb; ++j) {
                    omegas[j] = omegac * std::tan(phys::math::halfpi * ((j + 1.0f) / (Dimension::Nb + 1.0f)));
                    coeffs[j] = sqrt(2 * lambda / (Dimension::Nb + 1.0f)) * omegas[j];
                }
                break;
            }
            case HarmonicBathPolicy::Ohmic: {
                for (int j = 0; j < Dimension::Nb; ++j) {
                    omegas[j] = -omegac * std::log(1.0f - (j + 1.0f) / (Dimension::Nb + 1.0f));
                    coeffs[j] = sqrt(2 * lambda / (Dimension::Nb + 1.0f)) * omegas[j];
                }
                break;
            }
            case HarmonicBathPolicy::Closure: {
                // generate for eHEOM
                break;
            }
            case HarmonicBathPolicy::ReadFormula: {
                break;
            }
            case HarmonicBathPolicy::ReadFile:  // #b850, #pbi, #rub
            default: {
                try {
                    std::string bath_file =
                        _param->get_string({"model.bath_file", "model.bath.file"}, LOC(), "bath.dat");
                    std::ifstream ifs(bath_file);

                    std::string firstline, unit_str1, unit_str2, DISC_FLAG;
                    getline(ifs, firstline);
                    std::stringstream sstr(firstline);
                    sstr >> unit_str1 >> unit_str2;

                    phys::uval uw = phys::us::parse(unit_str1);
                    assert(uw.dim == phys::energy_d);
                    double w_unit = phys::us::conv(phys::au::unit, uw);

                    phys::uval uc = phys::us::parse(unit_str2);
                    assert(uc.dim == phys::energy_d);
                    double c_unit = phys::us::conv(phys::au::unit, uc);

                    // read data lines
                    double val;
                    for (int j = 0; j < Dimension::Nb; ++j) {
                        ifs >> DISC_FLAG;
                        if (DISC_FLAG == "WC") {                                   // input coefficients
                            if (ifs >> val) omegas[j] = val / w_unit;              ///< omegac ~ [energy_d]
                            if (ifs >> val) coeffs[j] = val / pow(c_unit, 1.5e0);  ///< coeffs ~ [energy_d]**1.5
                        } else if (DISC_FLAG == "WS") {                            // input Huang-Rhys factor
                            if (ifs >> val) omegas[j] = val / w_unit;              ///< omegac ~ [energy_d]
                            if (ifs >> val) coeffs[j] = 2.0f * sqrt(0.5e0 * omegas[j] * omegas[j] * omegas[j] * val);
                        } else if (DISC_FLAG == "WG") {                // input g coupling factor
                            if (ifs >> val) omegas[j] = val / w_unit;  ///< omegac ~ [energy_d]
                            if (ifs >> val) coeffs[j] = 2.0f * sqrt(0.5e0 * omegas[j] * omegas[j] * omegas[j]) * val;
                        } else {
                            throw kids_error("unknown discretization scheme");
                        }
                    }
                } catch (std::runtime_error& e) { throw kids_error("read bath.dat fails"); }
            }
        }
        if (is_et_transform) {
            for (int i = Dimension::Nb - 1; i > 0; --i) {
                omegas[i] = omegas[i - 1];
                coeffs[i] = coeffs[i - 1];
            }
            double omega0  = _param->get_real({"model.omega0"}, LOC(), 3.5e-4);
            double lambda0 = _param->get_real({"model.lambda0"}, LOC(), 2.39e-2);
            double coeff0  = _param->get_real({"model.coeff0"}, LOC(), std::sqrt(0.5 * lambda0) * omega0);
            // reorganization dressing
            double w2 = omega0 * omega0;
            for (int j = 1; j < Nb; ++j) { w2 += (coeffs[j] * coeffs[j]) / (omegas[j] * omegas[j]); }
            omegas[0] = std::sqrt(w2);
            coeffs[0] = coeff0;
            for (int j = 0; j < Dimension::N; ++j) w[j] = omegas[j % Nb];

            for (int j = 0, jk = 0; j < Dimension::N; ++j) {
                for (int k = 0; k < Dimension::N; ++k, ++jk) {
                    Kmat[jk] = (j == k) ? w[j] * w[j] : 0.0e0;
                    if (j == 0 && k > 0) Kmat[jk] = coeffs[k];
                    if (j > 0 && k == 0) Kmat[jk] = coeffs[j];
                    Tmod[jk] = (j == k) ? 1.0e0 : 0.0e0;
                }
            }
            // note in this scheme: ww, K, Tmod is not EigenSolve problem
        } else {
            for (int j = 0; j < Dimension::N; ++j) w[j] = omegas[j % Nb];
            for (int j = 0, jk = 0; j < Dimension::N; ++j) {
                for (int k = 0; k < Dimension::N; ++k, ++jk) {
                    Kmat[jk] = (j == k) ? w[j] * w[j] : 0.0e0;
                    Tmod[jk] = (j == k) ? 1.0e0 : 0.0e0;
                }
            }
        }
    }

    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);
    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    for (int j = 0; j < Dimension::N; ++j) {
        /* note:
            for finite temperature: Qoverbeta = 0.5*freq / dtanh(0.5*beta*freq)
            for zero temperature:   Qoverbeta = 0.5*freq
         */
        double Qoverbeta = (is_classical) ? ((beta > 0) ? 1.0f / beta : 0.0f)
                                          : (0.5e0 * w[j] / (beta > 0 ? std::tanh(0.5e0 * beta * w[j]) : 1.0f));
        x_sigma[j]       = std::sqrt(Qoverbeta / (w[j] * w[j]));
        p_sigma[j]       = std::sqrt(Qoverbeta);
    }
    if (is_et_transform) {
        for (int j = 0; j < Dimension::N; ++j) {  //
            x0[j] = (j % Nb == 0) ? -coeffs[0] / (omegas[0] * omegas[0]) : 0.0e0;
        }
    }
}

};  // namespace PROJECT_NS

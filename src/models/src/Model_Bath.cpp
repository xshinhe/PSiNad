#include "kids/Model_Bath.h"

#include <algorithm>

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Model_Bath::getName() { return "Model_Bath"; }

int Model_Bath::getType() const { return utils::hash(FUNCTION_NAME); }

double Model_Bath::J_Debye(double w) { return 2 * lambda * omegac * w / (w * w + omegac * omegac); }

double Model_Bath::J_Ohmic(double w) { return phys::math::pi * lambda / omegac * w * std::exp(-w / omegac); }

double Model_Bath::J(double w, double* w_arr, double* c_arr, int Nb) {
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

int Model_Bath::fun_Cw(kids_complex* Cw, double* w, int Nw, double* w_arr, double* c_arr, double beta, int Nb) {
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

void Model_Bath::setInputParam_impl(std::shared_ptr<Param> PM) {
    // size information
    Nb = _param->get_int({"model.Nb", "model.bath.Nb"}, LOC());
    assert(Nb == Dimension::Nb);

    bath_type         = BathPolicy::_from(_param->get_string({"model.bath_flag", "model.bath.flag"}, LOC(), "Debye"));
    omegac            = _param->get_real({"model.bath_omegac", "model.bath.omegac"}, LOC(), phys::energy_d, 1.0f);
    strength_type     = StrengthPolicy::_from(_param->get_string({"model.strength_flag"}, LOC(), "Lambda"));
    bath_is_classical = _param->get_bool({"model.bath_classic", "model.bath.classic"}, LOC(), false);

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

void Model_Bath::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    omegas = DS->def(DATA::model::bath::omegas);
    coeffs = DS->def(DATA::model::bath::coeffs);
    switch (bath_type) {
        case BathPolicy::Debye: {
            for (int j = 0; j < Dimension::Nb; ++j) {
                omegas[j] = omegac * std::tan(phys::math::halfpi * ((j + 1.0f) / (Dimension::Nb + 1.0f)));
                coeffs[j] = sqrt(2 * lambda / (Dimension::Nb + 1.0f)) * omegas[j];
            }
            break;
        }
        case BathPolicy::Ohmic: {
            for (int j = 0; j < Dimension::Nb; ++j) {
                omegas[j] = -omegac * std::log(1.0f - (j + 1.0f) / (Dimension::Nb + 1.0f));
                coeffs[j] = sqrt(2 * lambda / (Dimension::Nb + 1.0f)) * omegas[j];
            }
            break;
        }
        case BathPolicy::Closure: {
            break;
        }
        case BathPolicy::ReadFormula: {
            break;
        }
        case BathPolicy::ReadFile:  // #b850, #pbi, #rub
        default: {
            try {
                std::string   bath_file = _param->get_string({"model.bath_file", "model.bath.file"}, LOC(), "bath.dat");
                std::ifstream ifs(bath_file);

                std::string firstline, unit_str1, unit_str2, DIS_FLAG;
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
                    ifs >> DIS_FLAG;
                    if (DIS_FLAG == "WC") {                                    // input coefficients
                        if (ifs >> val) omegas[j] = val / w_unit;              ///< omegac ~ [energy_d]
                        if (ifs >> val) coeffs[j] = val / pow(c_unit, 1.5e0);  ///< coeffs ~ [energy_d]**1.5
                    } else if (DIS_FLAG == "WS") {                             // input Huang-Rhys factor
                        if (ifs >> val) omegas[j] = val / w_unit;              ///< omegac ~ [energy_d]
                        if (ifs >> val) coeffs[j] = 2.0f * sqrt(0.5e0 * omegas[j] * omegas[j] * omegas[j] * val);
                    } else if (DIS_FLAG == "WG") {                 // input g coupling factor
                        if (ifs >> val) omegas[j] = val / w_unit;  ///< omegac ~ [energy_d]
                        if (ifs >> val) coeffs[j] = 2.0f * sqrt(0.5e0 * omegas[j] * omegas[j] * omegas[j]) * val;
                    }
                }
            } catch (std::runtime_error& e) { throw kids_error("parse bath file error"); }
        }
    }

    x_sigma = DS->def(DATA::model::bath::x_sigma);
    p_sigma = DS->def(DATA::model::bath::p_sigma);
    for (int j = 0; j < Dimension::Nb; ++j) {
        /* note:
            for finite temperature: Qoverbeta = 0.5*freq / dtanh(0.5*beta*freq)
            for zero temperature:   Qoverbeta = 0.5*freq
         */
        double Qoverbeta = (bath_is_classical)
                               ? ((beta > 0) ? 1.0f / beta : 0.0f)
                               : (0.5e0 * omegas[j] / (beta > 0 ? std::tanh(0.5e0 * beta * omegas[j]) : 1.0f));
        x_sigma[j]       = std::sqrt(Qoverbeta / (omegas[j] * omegas[j]));
        p_sigma[j]       = std::sqrt(Qoverbeta);
    }
}

};  // namespace PROJECT_NS

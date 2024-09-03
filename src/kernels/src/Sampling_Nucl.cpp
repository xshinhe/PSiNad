#include "kids/Sampling_Nucl.h"

#include "kids/Kernel_Elec_Utils.h"
#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Representation.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

const std::string Sampling_Nucl::getName() { return "Sampling_Nulc"; }

int Sampling_Nucl::getType() const { return utils::hash(FUNCTION_NAME); }

void Sampling_Nucl::setInputParam_impl(std::shared_ptr<Param> PM) {
    sampling_type = NuclearSamplingPolicy::_from(
        _param->get_string({"solver.sampling_nuc_flag", "solver.sampling.nuc_flag"}, LOC(), "Gaussian"));
    sampling_file      = _param->get_string({"solver.sampling_file", "solver.sampling.file"}, LOC(), "init");
    screen_hfreq_type  = _param->get_int({"solver.screen_hfreq_type"}, LOC(), 0);
    ignore_nma         = _param->get_string({"solver.ignore_nma"}, LOC(), ",");
    double temperature = _param->get_real({"model.bath_temperature", "model.bath.temperature",  //
                                           "model.temperature", "solver.temperature"},
                                          LOC(), phys::temperature_d, 1.0f);
    beta               = 1.0f / (phys::au::k * temperature);  // please don't forget k_Boltzman
}

void Sampling_Nucl::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    x       = DS->def(DATA::integrator::x);
    p       = DS->def(DATA::integrator::p);
    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);
    w       = DS->def(DATA::model::w);
    mass    = DS->def(DATA::model::mass);
    Tmod    = DS->def(DATA::model::Tmod);
}

Status& Sampling_Nucl::initializeKernel_impl(Status& stat) { return stat; }

Status& Sampling_Nucl::executeKernel_impl(Status& stat) {
    // if restart, we should get all initial values and current values!
    // not that: all values defined by Kernel_Elec_Functions will also be recoveed later.
    // not that: all values defined by Kernel_Recorder will also be recovered later.
    if (_param->get_bool({"restart"}, LOC(), false)) {  //
        std::string loadfile = _param->get_string({"load"}, LOC(), "NULL");
        if (loadfile == "NULL" || loadfile == "" || loadfile == "null") loadfile = "restart.ds";
        _dataset->def(DATA::init::x, loadfile);
        _dataset->def(DATA::init::p, loadfile);
        _dataset->def(DATA::integrator::x, loadfile);
        _dataset->def(DATA::integrator::p, loadfile);
        return stat;
    }

    for (int iP = 0; iP < Dimension::P; ++iP) {  /// use P instead of P_NOW
        auto x = this->x.subspan(iP * Dimension::N, Dimension::N);
        auto p = this->p.subspan(iP * Dimension::N, Dimension::N);

        switch (sampling_type) {
            case NuclearSamplingPolicy::Fix: {
                for (int j = 0; j < Dimension::N; ++j) x[j] = x0[j], p[j] = 0.0e0;
                break;
            }
            case NuclearSamplingPolicy::Fix2: {
                for (int j = 0; j < Dimension::N; ++j) x[j] = x0[j], p[j] = p0[j];
                break;
            }
            case NuclearSamplingPolicy::ClassicalHO:
            case NuclearSamplingPolicy::WignerHO:
            case NuclearSamplingPolicy::QcHO: {
                Kernel_Random::rand_gaussian(x.data(), Dimension::N);
                Kernel_Random::rand_gaussian(p.data(), Dimension::N);
                for (int j = 0; j < Dimension::N; ++j) {
                    double Qoverbeta = (sampling_type == NuclearSamplingPolicy::ClassicalHO)
                                           ? ((beta > 0) ? 1.0f / beta : 0.0f)
                                           : (0.5e0 * w[j] / (beta > 0 ? std::tanh(0.5e0 * beta * w[j]) : 1.0f));
                    x_sigma[j]       = std::sqrt(Qoverbeta / (w[j] * w[j]));
                    p_sigma[j]       = std::sqrt(Qoverbeta);
                    x[j]             = x[j] * x_sigma[j];
                    p[j]             = p[j] * p_sigma[j];
                    if (sampling_type == NuclearSamplingPolicy::QcHO) {
                        double scale = sqrt(w[j] / (p[j] * p[j] + w[j] * w[j] * x[j] * x[j]));
                        x[j]         = x[j] * scale;
                        p[j]         = p[j] * scale;
                    }
                }
                break;
            }
            case NuclearSamplingPolicy::ClassicalNMA:
            case NuclearSamplingPolicy::WignerNMA:
            case NuclearSamplingPolicy::QcNMA: {
                double E_sampled_NMA = 0.0e0;
                Kernel_Random::rand_gaussian(x.data(), Dimension::N);
                Kernel_Random::rand_gaussian(p.data(), Dimension::N);
                for (int j = 0; j < Dimension::N; ++j) {
                    std::string findflag = utils::concat(",", j, ",");
                    bool find_in_ignore = ignore_nma.find(findflag) != std::string::npos;
                    if (find_in_ignore) {
                        std::cout << "here ignore nma of " << j << "-th frequency\n";
                    }
                    if (w[j] <= 0.0 || std::fabs(w[j]) < 1.0e-8 || find_in_ignore) {
                        x[j] = 0.0f, p[j] = 0.0f;
                        continue;
                    }
                    double Hratio = 0.0e0;
                    for(int k = 0; k < Dimension::N; ++k) {
                        if(mass[k] < 2000.0) Hratio += Tmod[k * Dimension::N + j] * Tmod[k * Dimension::N + j];
                    }
                    double Qoverbeta = (sampling_type == NuclearSamplingPolicy::ClassicalNMA)
                                           ? ((beta > 0) ? 1.0f / beta : 0.0f)
                                           : (0.5e0 * w[j] / (beta > 0 ? std::tanh(0.5e0 * beta * w[j]) : 1.0f));

                    double wincm = w[j] * phys::au_2_wn;
                    if(Hratio > 0.8 && wincm > 2500 && screen_hfreq_type > 0){
                        std::cout << "WARNING: using with screen on H freq: type = " << screen_hfreq_type << "\n";
                        if(screen_hfreq_type == 1){
                            x[j] = 0.0e0, p[j] = 0.0e0;
                        }
                        if(screen_hfreq_type == 2){
                            x[j] *= 0.1, p[j] *= 0.1;
                        }
                        if(screen_hfreq_type == 3){
                            Qoverbeta =  ((beta > 0) ? 1.0f / beta : 0.0f);
                            x[j] *= 0.1, p[j] *= 0.1;
                        }
                    }
                    E_sampled_NMA += w[j] * (x[j]*x[j] + p[j]*p[j]);

                    x_sigma[j]       = std::sqrt(Qoverbeta / (w[j] * w[j]));
                    p_sigma[j]       = std::sqrt(Qoverbeta);
                    x[j] *= x_sigma[j];
                    p[j] *= p_sigma[j];
                    if (sampling_type == NuclearSamplingPolicy::QcNMA) {
                        double scale = sqrt(w[j] / (p[j] * p[j] + w[j] * w[j] * x[j] * x[j]));
                        x[j]         = x[j] * scale;
                        p[j]         = p[j] * scale;
                    }
                }
                // save E_sampled_NMA to dataset;

                // transfrom normal-mode to cartesian coordinates (ARRAY should be .eval()!)
                ARRAY_MATMUL(x.data(), Tmod.data(), x.data(), Dimension::N, Dimension::N, 1);
                ARRAY_MATMUL(p.data(), Tmod.data(), p.data(), Dimension::N, Dimension::N, 1);
                for (int j = 0; j < Dimension::N; ++j) {
                    x[j] = x[j] / std::sqrt(mass[j]) + x0[j];
                    p[j] = p[j] * std::sqrt(mass[j]) + p0[j];
                }
                break;
            }
            case NuclearSamplingPolicy::Gaussian: {
                // x0 & p0 & x_sigma & p_sigma are provided
                // PRINT_ARRAY(x0, 1, Dimension::N);
                // PRINT_ARRAY(p0, 1, Dimension::N);
                // PRINT_ARRAY(x_sigma, 1, Dimension::N);
                // PRINT_ARRAY(p_sigma, 1, Dimension::N);
                Kernel_Random::rand_gaussian(x.data(), Dimension::N);
                Kernel_Random::rand_gaussian(p.data(), Dimension::N);
                for (int j = 0; j < Dimension::N; ++j) {
                    x[j] = x0[j] + x[j] * x_sigma[j];
                    p[j] = p0[j] + p[j] * p_sigma[j];
                }
                break;
            }
            case NuclearSamplingPolicy::ReadDataSet: {
                std::string open_file = sampling_file;
                if (!isFileExists(sampling_file)) open_file = utils::concat(sampling_file, stat.icalc, ".ds");
                std::string   stmp, eachline;
                std::ifstream ifs(open_file);
                while (getline(ifs, eachline)) {
                    if (eachline.find("init.x") != eachline.npos) {
                        getline(ifs, eachline);
                        for (int i = 0; i < Dimension::N; ++i) ifs >> x[i];
                    }
                    if (eachline.find("init.p") != eachline.npos) {
                        getline(ifs, eachline);
                        for (int i = 0; i < Dimension::N; ++i) ifs >> p[i];
                    }
                }
		//PRINT_ARRAY(x, 1, Dimension::N);
		//PRINT_ARRAY(p, 1, Dimension::N);
                //exit(0);
                break;
            }
            case NuclearSamplingPolicy::ReadXYZ: {
                //
                break;
            }
        }
    }
    _dataset->def(DATA::init::x, x);
    _dataset->def(DATA::init::p, p);
    return stat;
}

};  // namespace PROJECT_NS

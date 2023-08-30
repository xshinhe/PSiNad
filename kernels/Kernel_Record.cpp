#include "Kernel_Record.h"

#include "../core/Formula.h"
#include "../core/linalg.h"

namespace PROJECT_NS {

void Result::save(const std::string& fname, int ibegin, int length, bool with_header) {
    if (size <= 0) return;
    if (with_header) {  // print header
        if (ofs.is_open()) ofs.close();
        ofs.open(fname);
        ofs << FMT(8) << "t";
        ofs << FMT(8) << "stat";
        for (auto& v : header) ofs << FMT(8) << v;
        ofs << std::endl;
    }

    // ofs must has been open in this case
    if (!ofs.is_open()) throw std::runtime_error("result ofs is not open!");

    ibegin   = std::max({ibegin, 0});
    int iend = (length < 0) ? frame : std::min({ibegin + length, frame});
    for (int iframe = ibegin, idata = ibegin * size; iframe < iend; ++iframe) {
        ofs << FMT(8) << t0 + iframe * dt;
        ofs << FMT(8) << stat[iframe];
        for (int i = 0; i < size; ++i, ++idata) { ofs << FMT(8) << data[idata]; }
        ofs << std::endl;
    }
    ofs.flush();
}

Result::~Result() {
    if (ofs.is_open()) ofs.close();
}

Kernel_Record::~Kernel_Record() {
    if (ofs_samp.is_open()) ofs_samp.close();
    if (ofs_corr.is_open()) ofs_corr.close();
}

void Kernel_Record::read_param_impl(Param* PM) {
    dt        = PM->get<double>("dt", LOC(), phys::time_d);
    t0        = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    time_unit = PM->get<double>("time_unit", LOC(), phys::time_d, 1.0f);
    trace     = PM->get<bool>("trace", LOC(), false);
    directory = PM->get<std::string>("directory", LOC(), "default");
}

void Kernel_Record::init_data_impl(DataSet* DS) {
    istep_ptr = DS->reg<int>("timer.istep");
    sstep_ptr = DS->reg<int>("timer.sstep");
    isamp_ptr = DS->reg<int>("timer.isamp");
    nsamp_ptr = DS->reg<int>("timer.nsamp");
}

void Kernel_Record::init_calc_impl(int stat) {
    bool not_parsed = (Sampling_ID.size() == 0 && Correlation_ID1.size() == 0);
    auto& json      = *(_Param->pjson());
    if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
        Result& sampling    = get_sampling();
        Result& correlation = get_correlation();
        sampling.size       = 0;
        sampling.frame      = 1;
        correlation.size    = 0;
        correlation.frame   = (*nsamp_ptr);

        for (auto& j : (json["result"])) {
            if (!j.is_array()) continue;
            switch (j.size()) {
                case 1: {  // 1 point sampling
                    std::string f_str = j[0].get<std::string>();
                    int id            = Formula::regis_Formula(f_str, _DataSet, "integrator");
                    Sampling_ID.push_back(id);
                    auto& f = Formula::GLOBAL[id];
                    for (int i = 0; i < f.get_size(); ++i) {
                        if (f.get_res_type() == DataSet::Type::Real) {
                            sampling.header.push_back(utils::concat(f.name(), i));
                            sampling.size++;
                        } else {
                            sampling.header.push_back(utils::concat("R", f.name(), i));
                            sampling.header.push_back(utils::concat("I", f.name(), i));
                            sampling.size += 2;
                        }
                    }
                    break;
                }
                case 2: {  // 2 points correlation
                    std::string f_str1 = j[0].get<std::string>();
                    int id1            = Formula::regis_Formula(f_str1, _DataSet, "init");
                    Correlation_ID1.push_back(id1);

                    std::string f_str2 = j[1].get<std::string>();
                    int id2            = Formula::regis_Formula(f_str2, _DataSet, "integrator");
                    Correlation_ID2.push_back(id2);

                    auto& f1 = Formula::GLOBAL[id1];
                    auto& f2 = Formula::GLOBAL[id2];

                    for (int i1 = 0; i1 < f1.get_size(); ++i1) {
                        for (int i2 = 0; i2 < f2.get_size(); ++i2) {
                            if (f1.get_res_type() == DataSet::Type::Real && f2.get_res_type() == DataSet::Type::Real) {
                                correlation.header.push_back(utils::concat(f1.name(), i1, f2.name(), i2));
                                correlation.size++;
                            } else {
                                correlation.header.push_back(utils::concat("R", f1.name(), i1, f2.name(), i2));
                                correlation.header.push_back(utils::concat("I", f1.name(), i1, f2.name(), i2));
                                correlation.size += 2;
                            }
                        }
                    }
                    break;
                }
            }
        }
        sampling.t0 = t0;
        sampling.dt = dt * (*sstep_ptr) / time_unit;
        sampling.stat.resize(sampling.frame);
        sampling.data.resize(sampling.size * sampling.frame);
        correlation.t0 = t0;
        correlation.dt = dt * (*sstep_ptr) / time_unit;
        correlation.stat.resize(correlation.frame);
        correlation.data.resize(correlation.size * correlation.frame);
    }
}

int Kernel_Record::exec_kernel_impl(int stat) {
    Result& sampling      = get_sampling();
    Result& correlation   = get_correlation();
    bool do_record        = ((*istep_ptr) % (*sstep_ptr) == 0);
    bool do_record_header = ((*istep_ptr) == 0);

    if (do_record) {
        int sampling_idx0    = ((*isamp_ptr) % sampling.frame) * sampling.size;
        int correlation_idx0 = ((*isamp_ptr) % correlation.frame) * correlation.size;

        // calculate sampling
        for (int i = 0, idx = sampling_idx0; i < Sampling_ID.size(); ++i) {
            auto& f = Formula::GLOBAL[Sampling_ID[i]];
            for (int is = 0; is < f.get_size(); ++is) {
                switch (f.get_res_type()) {
                    case DataSet::Type::Real: {
                        sampling.data[idx++] = f.eval<num_real>(is);
                        break;
                    }
                    case DataSet::Type::Complex: {
                        num_complex v        = f.eval<num_complex>(is);
                        sampling.data[idx++] = std::real(v);
                        sampling.data[idx++] = std::imag(v);
                        break;
                    }
                }
            }
        }

        // calculate correlation
        for (int i = 0, idx = correlation_idx0; i < Correlation_ID1.size(); ++i) {
            auto& f1 = Formula::GLOBAL[Correlation_ID1[i]];
            auto& f2 = Formula::GLOBAL[Correlation_ID2[i]];
            for (int is1 = 0; is1 < f1.get_size(); ++is1) {
                for (int is2 = 0; is2 < f2.get_size(); ++is2) {
                    switch (f1.get_res_type()) {
                        case DataSet::Type::Real: {
                            num_real v1 = f1.eval<num_real>(is1);
                            switch (f2.get_res_type()) {
                                case DataSet::Type::Real: {
                                    num_real v2             = f2.eval<num_real>(is2);
                                    correlation.data[idx++] = v1 * v2;
                                    break;
                                }
                                case DataSet::Type::Complex: {
                                    num_complex v2          = f2.eval<num_complex>(is2);
                                    correlation.data[idx++] = std::real(v1 * v2);
                                    correlation.data[idx++] = std::imag(v1 * v2);
                                    break;
                                }
                            }
                            break;
                        }
                        case DataSet::Type::Complex: {
                            num_complex v1 = f1.eval<num_complex>(is1);
                            switch (f2.get_res_type()) {
                                case DataSet::Type::Real: {
                                    num_real v2             = f2.eval<num_real>(is2);
                                    correlation.data[idx++] = std::real(v1 * v2);
                                    correlation.data[idx++] = std::imag(v1 * v2);
                                    break;
                                }
                                case DataSet::Type::Complex: {
                                    num_complex v2          = f2.eval<num_complex>(is2);
                                    correlation.data[idx++] = std::real(v1 * v2);
                                    correlation.data[idx++] = std::imag(v1 * v2);
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }

        if (ARRAY_ISFINITE(sampling.data.data() + sampling_idx0, sampling.size) &&
            ARRAY_ISFINITE(correlation.data.data() + correlation_idx0, correlation.size)) {
            sampling.stat[(*isamp_ptr) % sampling.frame]       = 1;
            correlation.stat[(*isamp_ptr) % correlation.frame] = 1;
        } else {
            sampling.stat[(*isamp_ptr) % sampling.frame]       = 0;
            correlation.stat[(*isamp_ptr) % correlation.frame] = 0;
            Kernel::BREAK                                      = true;
        }

        if (trace) sampling.save(utils::concat(directory, "/samp", stat, ".dat"), 0, 1, do_record_header);
        // if (trace) correlation.save(utils::concat("corr", stat, ".dat"), *isamp_ptr, *isamp_ptr + 1,
        // do_record_header);
    }
    return 0;
}

Result Kernel_Record::_sampling;
Result Kernel_Record::_correlation;

};  // namespace PROJECT_NS

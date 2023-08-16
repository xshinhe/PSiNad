#include "Kernel_DataSetHandles.h"

namespace PROJECT_NS {

void Kernel_Load_DataSet::read_param_impl(Param* P) { fn = P->get<std::string>("load", LOC(), "NULL"); }

void Kernel_Load_DataSet::init_data_impl(DataSet* DS) { pDS = DS; }

int Kernel_Load_DataSet::exec_kernel_impl(int stat) {
    if (fn != "NULL") {
        try {
            pDS->load(fn);
        } catch (std::runtime_error& e) { throw state_load_error(fn); }
    }
    return 0;
}


void Kernel_Dump_DataSet::read_param_impl(Param* P) { fn = P->get<std::string>("dump", LOC(), "final"); }

void Kernel_Dump_DataSet::init_data_impl(DataSet* S) { pDS = S; }

int Kernel_Dump_DataSet::exec_kernel_impl(int stat) {
    if (fn == "null") return 0;
    try {
        pDS->dump(utils::concat(fn, stat, ".ds"));
    } catch (std::runtime_error& e) { throw state_dump_error(fn); }
    return 0;
}

void Result::save(const std::string& fname, double t0, double dt, bool with_header) {
    if (size <= 0) return;
    if (with_header) {  // print header
        if (ofs.is_open()) ofs.close();
        ofs.open(fname);
        ofs << FMT(8) << "t";
        for (auto& v : header) ofs << FMT(8) << v;
        ofs << std::endl;
    }
    // ofs must has been open in this case
    for (int iframe = 0, idata = 0; iframe < frame; ++iframe) {
        ofs << FMT(8) << t0 + iframe * dt;
        for (int i = 0; i < size; ++i, ++idata) { ofs << FMT(8) << data[idata]; }
        ofs << std::endl;
    }
}

Result::~Result() {
    if (ofs.is_open()) ofs.close();
}


Kernel_Record::~Kernel_Record() {
    if (ofs_samp.is_open()) ofs_samp.close();
    if (ofs_corr.is_open()) ofs_corr.close();
}

void Kernel_Record::read_param_impl(Param* P) {
    dt    = P->get<double>("dt", LOC(), phys::time_d);
    trace = P->get<bool>("trace", LOC(), false);
}

void Kernel_Record::init_data_impl(DataSet* DS) {
    istep_ptr = DS->reg<int>("timer.istep");
    sstep_ptr = DS->reg<int>("timer.sstep");
    isamp_ptr = DS->reg<int>("timer.isamp");
    nsamp_ptr = DS->reg<int>("timer.nsamp");

    Result& sampling    = get_sampling();
    Result& correlation = get_correlation();
    sampling.size       = 0;
    sampling.frame      = (trace) ? 1 : (*nsamp_ptr);
    correlation.size    = 0;
    correlation.frame   = (trace) ? 1 : (*nsamp_ptr);

    auto& json = *(_Param->pjson());
    if (json.count("result") == 1 && json["result"].is_array()) {
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
    }
    sampling.data.resize(sampling.size * sampling.frame);
    correlation.data.resize(correlation.size * correlation.frame);
}

void Kernel_Record::init_calc_impl(int stat){};
//  {
//     bool not_parsed = (Sampling_ID.size() == 0 && Correlation_ID1.size() == 0);
//     auto& json      = *(_Param->pjson());
//     if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
//         Result& sampling    = get_sampling();
//         Result& correlation = get_correlation();
//         sampling.size       = 0;
//         sampling.frame      = (trace) ? 1 : (*nsamp_ptr);
//         correlation.size    = 0;
//         correlation.frame   = (trace) ? 1 : (*nsamp_ptr);

//         for (auto& j : (json["result"])) {
//             if (!j.is_array()) continue;
//             switch (j.size()) {
//                 case 1: {  // 1 point sampling
//                     std::string f_str = j[0].get<std::string>();
//                     int id            = Formula::regis_Formula(f_str, _DataSet, "integrator");
//                     Sampling_ID.push_back(id);
//                     auto& f = Formula::GLOBAL[id];
//                     for (int i = 0; i < f.get_size(); ++i) {
//                         if (f.get_res_type() == DataSet::Type::Real) {
//                             sampling.header.push_back(utils::concat(f.name(), i));
//                             sampling.size++;
//                         } else {
//                             sampling.header.push_back(utils::concat("R", f.name(), i));
//                             sampling.header.push_back(utils::concat("I", f.name(), i));
//                             sampling.size += 2;
//                         }
//                     }
//                     break;
//                 }
//                 case 2: {  // 2 points correlation
//                     std::string f_str1 = j[0].get<std::string>();
//                     int id1            = Formula::regis_Formula(f_str1, _DataSet, "init");
//                     Correlation_ID1.push_back(id1);

//                     std::string f_str2 = j[1].get<std::string>();
//                     int id2            = Formula::regis_Formula(f_str2, _DataSet, "integrator");
//                     Correlation_ID2.push_back(id2);

//                     auto& f1 = Formula::GLOBAL[id1];
//                     auto& f2 = Formula::GLOBAL[id2];

//                     for (int i1 = 0; i1 < f1.get_size(); ++i1) {
//                         for (int i2 = 0; i2 < f2.get_size(); ++i2) {
//                             if (f1.get_res_type() == DataSet::Type::Real && f2.get_res_type() == DataSet::Type::Real)
//                             {
//                                 correlation.header.push_back(utils::concat(f1.name(), i1, f2.name(), i2));
//                                 correlation.size++;
//                             } else {
//                                 correlation.header.push_back(utils::concat("R", f1.name(), i1, f2.name(), i2));
//                                 correlation.header.push_back(utils::concat("I", f1.name(), i1, f2.name(), i2));
//                                 correlation.size += 2;
//                             }
//                         }
//                     }
//                     break;
//                 }
//             }
//         }
//         sampling.data.resize(sampling.size * sampling.frame);
//         correlation.data.resize(correlation.size * correlation.frame);
//     }
// }

int Kernel_Record::exec_kernel_impl(int stat) {
    Result& sampling      = get_sampling();
    Result& correlation   = get_correlation();
    bool do_record        = ((*istep_ptr) % (*sstep_ptr) == 0);
    bool do_record_header = ((*istep_ptr) == 0);

    if (do_record) {
        double time          = (*istep_ptr) * dt;
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

        if (trace) sampling.save(utils::concat("samp", stat, ".dat"), time, 0, do_record_header);

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

        if (trace) correlation.save(utils::concat("corr", stat, ".dat"), time, 0, do_record_header);
    }
    return 0;
}

Result Kernel_Record::_sampling;
Result Kernel_Record::_correlation;

};  // namespace PROJECT_NS

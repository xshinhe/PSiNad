#include "Kernel_Record.h"

#include "../core/Formula.h"
#include "../core/linalg.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(8) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

namespace PROJECT_NS {

void Result::stack(kids_dtype itype, std::string str) {
    switch (itype) {
        case kids_real_type:
            header.push_back(str);
            size++;
            break;
        case kids_complex_type:
            header.push_back(utils::concat("R", str));
            header.push_back(utils::concat("I", str));
            size += 2;
            break;
    }
}

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

    if (false && with_header) {  // print header
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

    if (false) {  // add filter to select which data should be outputed here! @TODO
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
}

Result::~Result() {
    if (ofs.is_open()) ofs.close();
}

Kernel_Record::~Kernel_Record() {
    if (ofs_corr.is_open()) ofs_corr.close();
}

void Kernel_Record::read_param_impl(Param* PM) {
    dt        = PM->get<double>("dt", LOC(), phys::time_d);
    t0        = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    time_unit = PM->get<double>("time_unit", LOC(), phys::time_d, 1.0f);
    directory = PM->get<std::string>("directory", LOC(), "default");
}

void Kernel_Record::init_data_impl(DataSet* DS) {
    istep_ptr                     = DS->def<kids_int>("iter.istep");
    sstep_ptr                     = DS->def<kids_int>("iter.sstep");
    isamp_ptr                     = DS->def<kids_int>("iter.isamp");
    nsamp_ptr                     = DS->def<kids_int>("iter.nsamp");
    at_samplingstep_initially_ptr = DS->def<kids_bool>("iter.at_samplingstep_initially");
}

inline std::string contacted_hdr(const std::string& s1, int i1, const std::string& s2, int i2) {
    return utils::concat(s1, i1, s2, i2);
}

inline std::string hyphen_hdr(const std::string& s1, int i1, const std::string& s2, int i2) {
    return utils::concat(s1, i1, "-", s2, i2);
}

inline std::string bracket_hdr(const std::string& s1, int i1, const std::string& s2, int i2) {
    return utils::concat(s1, s2, "[", i1, ",", i2, "]");
}

void Kernel_Record::token_array(Result& correlation, Param::JSON& j) {  // @deprecated[2024] old format
    Record_Item item_tmp;

    switch (j.size()) {
        case 1: {  // 1 point sampling
            item_tmp.v0      = "1";
            item_tmp.vt      = j[0].get<std::string>();
            item_tmp.fml_ID0 = Formula::regis_Formula(item_tmp.v0, _DataSet, "init");
            item_tmp.fml_IDt = Formula::regis_Formula(item_tmp.vt, _DataSet, "integrator");
            item_tmp.name    = "_";
            item_tmp.save    = "samp";
            // item_tmp.ofs_ID  = ofsm.push(item_tmp.save);

            Record_List.push_back(item_tmp);

            auto& f = Formula::GLOBAL[item_tmp.fml_IDt];
            for (int i = 0; i < f.get_size(); ++i) {  //
                correlation.stack(f.get_res_type(), utils::concat(f.name(), i));
            }
            break;
        }
        case 2: {  // 2 points correlation
            item_tmp.v0      = j[0].get<std::string>();
            item_tmp.vt      = j[1].get<std::string>();
            item_tmp.fml_ID0 = Formula::regis_Formula(item_tmp.v0, _DataSet, "init");
            item_tmp.fml_IDt = Formula::regis_Formula(item_tmp.vt, _DataSet, "integrator");
            item_tmp.name    = "_";
            item_tmp.save    = "corr";
            // item_tmp.ofs_ID  = ofsm.push(item_tmp.save);

            Record_List.push_back(item_tmp);

            auto&      f1     = Formula::GLOBAL[item_tmp.fml_ID0];
            auto&      f2     = Formula::GLOBAL[item_tmp.fml_IDt];
            kids_dtype f_type = (f1.get_res_type() == kids_real_type && f2.get_res_type() == kids_real_type)
                                    ? kids_real_type
                                    : kids_complex_type;

            /**
             * Output format for correlation: A -- B -- C -- D correlation
             * 1) A1B2C3D4: contacted format
             * 1) A1-B2-C3-D4: hyphen format
             * 2) ABCD[1,2,3,4]: bracket format
             */
            for (int i1 = 0; i1 < f1.get_size(); ++i1) {
                for (int i2 = 0; i2 < f2.get_size(); ++i2) {
                    correlation.stack(f_type, contacted_hdr(f1.name(), i1, f2.name(), i2));
                }
            }
            break;
        }
    }
}

void Kernel_Record::token_object(Result& correlation, Param::JSON& j) {  // @recomment new format
    Record_Item item_tmp;
    item_tmp.v0      = (j.count("v0") == 1) ? j["v0"].get<std::string>() : "1";
    item_tmp.vt      = (j.count("vt") == 1) ? j["vt"].get<std::string>() : "1";
    item_tmp.name    = (j.count("name") == 1) ? j["name"].get<std::string>() : "_";
    item_tmp.save    = (j.count("save") == 1) ? j["save"].get<std::string>() : "corr";
    item_tmp.fml_ID0 = Formula::regis_Formula(item_tmp.v0, _DataSet, "init");
    item_tmp.fml_IDt = Formula::regis_Formula(item_tmp.vt, _DataSet, "integrator");
    item_tmp.ofs_ID  = ofsm.push(item_tmp.save);

    Record_List.push_back(item_tmp);

    auto& f1 = Formula::GLOBAL[item_tmp.fml_ID0];
    auto& f2 = Formula::GLOBAL[item_tmp.fml_IDt];

    kids_dtype f_type = (f1.get_res_type() == kids_real_type && f2.get_res_type() == kids_real_type)
                            ? kids_real_type
                            : kids_complex_type;
    for (int i1 = 0; i1 < f1.get_size(); ++i1) {
        for (int i2 = 0; i2 < f2.get_size(); ++i2) {
            correlation.stack(f_type, utils::concat(f1.name(), i1, f2.name(), i2));
        }
    }
}

void Kernel_Record::init_calc_impl(int stat) {
    bool  not_parsed = Record_List.size() == 0;
    auto& json       = *(_Param->pjson());
    if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
        Result& correlation = get_correlation();
        correlation.size    = 0;
        correlation.frame   = (*nsamp_ptr);

        // add default record here
        // Param::JSON def_j = Param::JSON::parse("[\"K0\", \"K0\"]");
        // token_array(correlation, def_j);

        for (auto& j : (json["result"])) {
            if (j.is_array()) token_array(correlation, j);
            if (j.is_object()) token_object(correlation, j);
        }
        correlation.t0 = t0;
        correlation.dt = sstep_ptr[0] * dt / time_unit;
        correlation.stat.resize(correlation.frame);
        correlation.data.resize(correlation.size * correlation.frame);
    }
}

int Kernel_Record::exec_kernel_impl(int stat) {
    Result& correlation = get_correlation();
    if (at_samplingstep_initially_ptr[0]) {
        int correlation_idx0 = ((*isamp_ptr) % correlation.frame) * correlation.size;
        for (int i = 0, idx = correlation_idx0; i < Record_List.size(); ++i) {
            auto& f1 = Formula::GLOBAL[Record_List[i].fml_ID0];
            auto& f2 = Formula::GLOBAL[Record_List[i].fml_IDt];

            for (int is1 = 0; is1 < f1.get_size(); ++is1) {
                for (int is2 = 0; is2 < f2.get_size(); ++is2) {
                    switch (f1.get_res_type()) {
                        case kids_real_type: {
                            kids_real v1 = f1.eval<kids_real>(is1);
                            switch (f2.get_res_type()) {
                                case kids_real_type: {
                                    kids_real v2            = f2.eval<kids_real>(is2);
                                    correlation.data[idx++] = v1 * v2;
                                    break;
                                }
                                case kids_complex_type: {
                                    kids_complex v2         = f2.eval<kids_complex>(is2);
                                    correlation.data[idx++] = std::real(v1 * v2);
                                    correlation.data[idx++] = std::imag(v1 * v2);
                                    break;
                                }
                            }
                            break;
                        }
                        case kids_complex_type: {
                            kids_complex v1 = f1.eval<kids_complex>(is1);
                            switch (f2.get_res_type()) {
                                case kids_real_type: {
                                    kids_real v2            = f2.eval<kids_real>(is2);
                                    correlation.data[idx++] = std::real(v1 * v2);
                                    correlation.data[idx++] = std::imag(v1 * v2);
                                    break;
                                }
                                case kids_complex_type: {
                                    kids_complex v2         = f2.eval<kids_complex>(is2);
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

        if (ARRAY_ISFINITE(correlation.data.data() + correlation_idx0, correlation.size)) {
            correlation.stat[(*isamp_ptr) % correlation.frame] = 1;
        } else {
            correlation.stat[(*isamp_ptr) % correlation.frame] = 0;
            Kernel::BREAK                                      = true;
        }
    }
    return 0;
}

Result Kernel_Record::_correlation;

};  // namespace PROJECT_NS

#include "Kernel_Record.h"

#include "../core/Einsum.h"
#include "../core/Formula.h"
#include "../core/linalg.h"
#include "../core/vars_list.h"

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
    istep_ptr                     = DS->def(DATA::iter::istep);
    sstep_ptr                     = DS->def(DATA::iter::sstep);
    isamp_ptr                     = DS->def(DATA::iter::isamp);
    nsamp_ptr                     = DS->def(DATA::iter::nsamp);
    at_samplingstep_initially_ptr = DS->def(DATA::iter::at_samplingstep_initially);
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

struct VarItem {
   public:
    VarItem(const std::string& input, DataSet* DS, bool def_in = true) : input{input} {
        std::string            str = input;
        std::string::size_type ipos;

        ipos = str.find(":");
        type = (ipos == std::string::npos) ? "" : str.substr(ipos + 1, str.size());
        str  = str.substr(0, ipos);

        ipos  = str.find("<");
        index = (ipos == std::string::npos) ? "" : str.substr(ipos + 1, str.size());
        index.erase(index.find_last_not_of(">") + 1);
        str = str.substr(0, ipos);

        ipos  = str.find("{");
        field = (ipos == std::string::npos) ? "" : str.substr(ipos + 1, str.size());
        field.erase(field.find_last_not_of("}") + 1);
        name = str.substr(0, ipos);

        ipos  = field.find("@");
        time  = (ipos == std::string::npos) ? "" : field.substr(ipos + 1, str.size());
        field = (ipos == std::string::npos) ? field : field.substr(0, ipos);

        if (field == "I") field = "integrator";
        if (field == "M") field = "model";
        if (field == "R") field = "result";
        if (field == "" && def_in) field = "integrator";
        if (field == "" && !def_in) field = "result";
        if (time == "" || time == "T") time = "t";
        if (def_in) def(DS);
    }

    void def(DataSet* DS) {
        std::string key                = utils::concat(field, ".", name);
        std::tie(_type, _data, _shape) = DS->obtain(key);
    }

    std::string input;
    std::string name;
    std::string field;
    std::string index;
    std::string type;
    std::string time;
    kids_dtype  _type;
    void*       _data;
    Shape*      _shape;
};

struct Record_Rule {
    Record_Rule(const std::string& rule, DataSet* DS, const std::string& mode = "average",
                const std::string& save = "res.dat")
        : rule{rule}, mode{mode}, save{save} {
        std::string res0_str;
        std::string vars_str;
        std::string item_str;

        // here "=" is the separator
        auto ipos = rule.find("=");
        vars_str  = rule.substr(0, ipos);
        expr_str  = rule.substr(ipos + 1, rule.size());
        expr_str.erase(0, expr_str.find_first_not_of(" "));
        expr_str.erase(expr_str.find_last_not_of(" ") + 1);

        // here "(,)" are the separators
        ipos     = vars_str.find("(");
        res0_str = vars_str.substr(0, ipos);
        res0_str.erase(0, res0_str.find_first_not_of(" "));
        res0_str.erase(res0_str.find_last_not_of(" ") + 1);
        _res.push_back(VarItem(res0_str, DS, false));

        vars_str = vars_str.substr(ipos + 1, vars_str.size());
        while ((ipos = vars_str.find(",")) != std::string::npos) {
            item_str = vars_str.substr(0, ipos);
            item_str.erase(0, item_str.find_first_not_of(" "));
            item_str.erase(item_str.find_last_not_of(" ") + 1);
            vars_str = vars_str.substr(ipos + 1, vars_str.size());
            vars.push_back(VarItem(item_str, DS));
        }
        ipos     = vars_str.find(")");
        item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars.push_back(VarItem(item_str, DS));

        std::vector<std::vector<std::size_t>> input_shapes;
        eins_str = "";
        for (int i = 0; i < vars.size(); ++i) {
            if (i != 0) eins_str += ",";
            eins_str += vars[i].index;
            input_shapes.push_back(vars[i]._shape->dims());
        }
        if (_res[0].index != "") {
            eins_str += "->";
            eins_str += _res[0].index;
        }
        for (auto&& i : input_shapes) {
            for (auto&& j : i) { std::cout << j << ","; }
            std::cout << "\n";
        }
        std::cout << eins_str << "\n";
        auto&& EH = EinsumHelper(eins_str, input_shapes);
        for (auto&& i : EH.dh_output.dims) {
            { std::cout << i << ","; }
        }
        std::cout << "\n";

        DS->def<kids_complex>(utils::concat(_res[0].field, ".", _res[0].name),  //
                              EH.dh_output.dims, "CUSTOM");
        _res[0].def(DS);
    }

    std::string          rule;
    std::string          mode;
    std::string          save;
    std::string          expr_str;
    std::string          eins_str;
    std::vector<VarItem> vars;
    std::vector<int>     fids;
    std::vector<VarItem> _res;
};

void Kernel_Record::token(Result& correlation, Param::JSON& j) {
    std::string rule, mode, save;
    if (j.is_string()) {
        rule = j.get<std::string>();
        mode = "average";
        save = "res.dat";
    } else if (j.is_object()) {
        rule = j["rule"].get<std::string>();
        mode = j["mode"].get<std::string>();
        save = j["save"].get<std::string>();
    } else {
        throw std::runtime_error("parse error");
    }
    Record_Rule rule_item(rule, _DataSet, mode, save);

    std::cout << rule_item.rule << "\n";
    std::cout << rule_item.mode << "\n";
    std::cout << rule_item.save << "\n";
    std::cout << rule_item.expr_str << "\n";
    std::cout << rule_item.eins_str << "\n";
    for (auto&& v : rule_item.vars) {
        std::cout << v.input << "\n";
        std::cout << v.name << "\n";
        std::cout << v.field << "\n";
        std::cout << v.index << "\n";
        std::cout << v.type << "\n";
        std::cout << v.time << "\n";
        std::cout << (int) v._type << "\n";
        for (auto&& id : (*(v._shape)).dims()) { std::cout << id << ","; }
        std::cout << "\n";
    }

    for (auto&& v : rule_item._res) {
        std::cout << v.input << "\n";
        std::cout << v.name << "\n";
        std::cout << v.field << "\n";
        std::cout << v.index << "\n";
        std::cout << v.type << "\n";
        std::cout << v.time << "\n";
        std::cout << (int) v._type << "\n";
        for (auto&& id : (*(v._shape)).dims()) { std::cout << id << ","; }
        std::cout << "\n";
    }
    exit(0);
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
            if (j.is_string()) token(correlation, j);
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

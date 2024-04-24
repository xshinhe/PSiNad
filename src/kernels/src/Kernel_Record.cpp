#include "kids/Kernel_Record.h"

#include "kids/Einsum.h"
#include "kids/Formula.h"
#include "kids/linalg.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

/** custumized einsum-like function operation
 * @tparam          T               data type
 * @param[in]       EH              EinsumHelper object
 * @param[in]       data_inputs     vector of pointers of data of input tensors
 * @param[inout]    data_output     pointer stored data of output tensor
 */
template <typename T, typename Tout>
void einsum_fun(EinsumHelper&                  EH,           //
                FPARSER<T>&                    FP,           //
                const std::vector<void*>&      data_inputs,  //
                const std::vector<kids_dtype>& data_idtype,  //
                Tout*                          data_output   //
) {
    auto& einsum_dims   = EH.einsum_dims;
    auto& einsum_iposes = EH.einsum_iposes;
    auto& ipos_inputs   = EH.ipos_inputs;
    // ipos_output
    auto& dh_inputs          = EH.dh_inputs;
    auto& dh_output_mapldims = EH.dh_output.mapldims;

    std::size_t total_loop   = EH.total_loop;
    std::size_t total_tensor = EH.total_tensor;
    std::size_t total_esidx  = EH.total_esidx;
    std::size_t imax         = EH.count3 - 1;
    std::size_t imin         = EH.count1;

    std::vector<T> holder(total_tensor);  // temporary vector

    memset(einsum_iposes.data(), 0, total_esidx * sizeof(std::size_t));
    memset(ipos_inputs.data(), 0, total_tensor * sizeof(std::size_t));
    data_output[0] = Tout(0);
    for (std::size_t iloop = 0, ipos_output = 0; iloop < total_loop; ++iloop) {
        for (int iten = 0; iten < total_tensor; ++iten) {  //
            switch (data_idtype[iten]) {
                case kids_int_type:
                    holder[iten] = cast_at<T, kids_int>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case kids_real_type:
                    holder[iten] = cast_at<T, kids_real>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case kids_complex_type:
                    holder[iten] = cast_at<T, kids_complex>(data_inputs[iten], ipos_inputs[iten]);
                    break;
            }
        }
        data_output[ipos_output] += cast<Tout>(FP.eval(holder.data()));

        std::size_t i = imax;
        while (++einsum_iposes[i] == einsum_dims[i] && i > imin) { einsum_iposes[i--] = 0; }

        for (int iten = 0; iten < total_tensor; ++iten)  //
            ipos_inputs[iten] += dh_inputs[iten].mapldims[i];

        ipos_output += dh_output_mapldims[i];
        if (i < EH.count2) { data_output[ipos_output] = Tout(0); }
    }
}

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

void Kernel_Record::setInputParam_impl(std::shared_ptr<Param>& PM) {
    dt        = PM->get_double("dt", LOC(), phys::time_d);
    t0        = PM->get_double("t0", LOC(), phys::time_d, 0.0f);
    time_unit = PM->get_double("time_unit", LOC(), phys::time_d, 1.0f);
    directory = PM->get_string("directory", LOC(), "default");
}

void Kernel_Record::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
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

        if (field == "0") field = "init";
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
                const std::string& save = "res.dat", int nsamp = 1)
        : rule{rule}, mode{mode}, save{save}, nsamp{nsamp} {
        std::string res0_str;
        std::string vars_str;
        std::string item_str;

        // here "=" is the separator
        auto ipos = rule.find("=");
        vars_str  = (ipos == std::string::npos) ? rule : rule.substr(0, ipos);
        expr_str  = (ipos == std::string::npos) ? "" : rule.substr(ipos + 1, rule.size());
        expr_str.erase(0, expr_str.find_first_not_of(" "));
        expr_str.erase(expr_str.find_last_not_of(" ") + 1);

        // here "(,)" are the separators
        ipos     = vars_str.find("(");
        res0_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(0, ipos);
        res0_str.erase(0, res0_str.find_first_not_of(" "));
        res0_str.erase(res0_str.find_last_not_of(" ") + 1);
        _res        = std::shared_ptr<VarItem>(new VarItem(res0_str, DS, false));
        _res->field = "result";

        vars_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(ipos + 1, vars_str.size());
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

        eins_str = "";
        for (int i = 0; i < vars.size(); ++i) {
            if (i != 0) eins_str += ",";
            eins_str += vars[i].index;
            input_shapes.push_back(vars[i]._shape->dims());
            data_inputs.push_back(vars[i]._data);
            data_idtype.push_back(vars[i]._type);
        }
        if (_res->index != "") eins_str += utils::concat("->", _res->index);
        EH = std::shared_ptr<EinsumHelper>(new EinsumHelper(eins_str, input_shapes));
        std::vector<std::size_t> shape_output_all;
        // dressed the shape_output_all
        shape_output_all.push_back(nsamp);
        for (auto&& i : EH->dh_output.dims) shape_output_all.push_back(i);
        //
        DS->def_complex(utils::concat(_res->field, ".", _res->name),  //
                        shape_output_all, "CUSTOM DEFINED");

        expr_type = kids_void_type;
        for (auto&& v : vars) {
            if (v._type != kids_int_type && v._type != kids_real_type && v._type != kids_complex_type)
                throw std::runtime_error("access bad type");
            if (v._type > expr_type) expr_type = v._type;
        }
        _res->def(DS);
        std::string varslist = "";
        for (int i = 0; i < vars.size(); ++i) {  //
            varslist += (i == 0) ? utils::concat("_", i) : utils::concat(",_", i);
        }
        if (expr_str == "" && vars.size() > 0) {
            for (int i = 0; i < vars.size(); ++i) {  //
                expr_str += (i == 0) ? utils::concat("_", i) : utils::concat(" * _", i);
            }
        }
        switch (expr_type) {
            case kids_real_type:
                expr_id = FPARSER<kids_real>::regis_FPARSER(expr_str, varslist);
                break;
            case kids_complex_type:
                expr_id = FPARSER<kids_complex>::regis_FPARSER(expr_str, varslist);
                break;
        }
    }

    void calc(int isamp = 0) {
        switch (expr_type) {
            case kids_real_type: {
                auto&& eval = (FPARSER<kids_real>::GLOBAL[expr_id]);
                if (_res->_type == kids_real_type) {
                    kids_real* resdata = (kids_real*) (_res->_data) + isamp * EH->total_loop;
                    einsum_fun(*EH, eval, data_inputs, data_idtype, resdata);
                }
                if (_res->_type == kids_complex_type) {
                    kids_complex* resdata = (kids_complex*) (_res->_data) + isamp * EH->total_loop;
                    einsum_fun(*EH, eval, data_inputs, data_idtype, resdata);
                }
                break;
            }
            case kids_complex_type: {
                auto&& eval = (FPARSER<kids_complex>::GLOBAL[expr_id]);
                if (_res->_type == kids_real_type) {
                    kids_real* resdata = (kids_real*) (_res->_data) + isamp * EH->total_loop;
                    einsum_fun(*EH, eval, data_inputs, data_idtype, resdata);
                }
                if (_res->_type == kids_complex_type) {
                    kids_complex* resdata = (kids_complex*) (_res->_data) + isamp * EH->total_loop;
                    einsum_fun(*EH, eval, data_inputs, data_idtype, resdata);
                }
                break;
            }
        }
    }

    int                                   nsamp;
    std::string                           rule;
    std::string                           mode;
    std::string                           save;
    std::string                           eins_str;
    std::shared_ptr<EinsumHelper>         EH;
    std::string                           expr_str;
    kids_dtype                            expr_type;
    std::size_t                           expr_id;
    std::vector<VarItem>                  vars;
    std::vector<std::vector<std::size_t>> input_shapes;
    std::vector<void*>                    data_inputs;
    std::vector<kids_dtype>               data_idtype;
    std::shared_ptr<VarItem>              _res;
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

    Rules.push_back(std::shared_ptr<Record_Rule>(new Record_Rule(rule, _DataSet, mode, save, nsamp_ptr[0])));

    auto&& rule_item = *(Rules[0]);

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

    std::cout << rule_item.expr_type << "\n";
    std::cout << rule_item._res->input << "\n";
    std::cout << rule_item._res->name << "\n";
    std::cout << rule_item._res->field << "\n";
    std::cout << rule_item._res->index << "\n";
    std::cout << rule_item._res->type << "\n";
    std::cout << rule_item._res->time << "\n";
    std::cout << (int) rule_item._res->_type << "\n";
    for (auto&& id : (*(rule_item._res->_shape)).dims()) { std::cout << id << ","; }
    std::cout << "\n\n";
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

Status& Kernel_Record::initializeKernel_impl(Status& stat) {
    bool  not_parsed = Record_List.size() == 0;
    auto& json       = *(_param->pjson());
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

Status& Kernel_Record::executeKernel_impl(Status& stat) {
    Result& correlation = get_correlation();
    if (at_samplingstep_initially_ptr[0]) {
        for (auto&& irecord : Rules) { irecord->calc(isamp_ptr[0]); }
    }
    if (at_samplingstep_initially_ptr[0] && false) {
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
            // Kernel::BREAK                                      = true;
        }
    }
    return 0;
}

Result Kernel_Record::_correlation;

};  // namespace PROJECT_NS

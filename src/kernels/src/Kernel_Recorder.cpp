#include "kids/Kernel_Recorder.h"

#include <algorithm>

#include "kids/Einsum.h"
#include "kids/RuleEvaluator.h"
#include "kids/RuleSet.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Recorder::getName() { return "Kernel_Recorder"; }

int Kernel_Recorder::getType() const { return utils::hash(FUNCTION_NAME); }

Kernel_Recorder::Kernel_Recorder() {
    _ruleset = std::shared_ptr<RuleSet>(new RuleSet());  //
}

Kernel_Recorder::~Kernel_Recorder(){};

void Kernel_Recorder::setInputParam_impl(std::shared_ptr<Param> PM) {
    dt        = _param->get_real({"model.dt", "solver.dt"}, LOC(), phys::time_d);
    t0        = _param->get_real({"model.t0", "solver.t0"}, LOC(), phys::time_d, 0.0f);
    time_unit = _param->get_real({"model.time_unit", "solver.time_unit"}, LOC(), phys::time_d, 1.0f);
    occ0      = _param->get_int({"model.occ", "solver.occ"}, LOC(), -1);
}

void Kernel_Recorder::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    istep_ptr                = DS->def(DATA::flowcontrol::istep);
    sstep_ptr                = DS->def(DATA::flowcontrol::sstep);
    isamp_ptr                = DS->def(DATA::flowcontrol::isamp);
    nsamp_ptr                = DS->def(DATA::flowcontrol::nsamp);
    kids_real* per_time_unit = DS->def(DATA::flowcontrol::pertimeunit);
    per_time_unit[0]         = 1.0e0 / time_unit;
}

void Kernel_Recorder::parse() {
    if (!_param->is_array("record")) return;

    int it = 0;
    while (true) {
        std::string firstkey = utils::concat("record.", it);
        if (!_param->has_key(firstkey)) break;

        std::string rule = "", mode = "average", save = "res.dat";

        if (false) {
        } else if (_param->is_object(firstkey)) {
            std::string mode_key = utils::concat(firstkey, ".", "mode");
            std::string save_key = utils::concat(firstkey, ".", "save");
            std::string rule_key = utils::concat(firstkey, ".", "rule");
            if (_param->has_key(mode_key)) mode = _param->get_string({mode_key}, LOC());
            if (_param->has_key(save_key)) save = _param->get_string({save_key}, LOC());
            if (_param->has_key(rule_key)) {
                rule = _param->get_string({rule_key}, LOC());
            } else {
                throw kids_error("parser rule error");
            }
        } else if (_param->is_array(firstkey)) {
            std::string firstv0key = utils::concat(firstkey, ".", "0");
            std::string firstv1key = utils::concat(firstkey, ".", "1");
            std::string v0, v1;
            if (_param->has_key(firstv0key)) v0 = _param->get_string({firstv0key}, LOC());
            if (_param->has_key(firstv1key)) v1 = _param->get_string({firstv1key}, LOC());

            auto ipos0 = v0.find_first_of("#");
            if (ipos0 != std::string::npos) v0 = v0.substr(0, ipos0);
            auto ipos1 = v1.find_first_of("#");
            if (ipos1 != std::string::npos) v1 = v1.substr(0, ipos1);
            rule = utils::concat("_", it, "(", v0, ",", v1, ")");
        } else if (_param->is_string(firstkey)) {
            rule = _param->get_string({firstkey}, LOC());
        } else {
            throw kids_error("unknown type");
        }

        if (std::find(opened_files.begin(), opened_files.end(), save) == opened_files.end()) {
            std::shared_ptr<RuleEvaluator> record_time_rule(  //
                new RuleEvaluator("time(t{flowcontrol}:R, pertimeunit{flowcontrol}:R)", _dataset, "copy", save,
                                  nsamp_ptr[0]));
            _ruleset->registerRules(record_time_rule);
            opened_files.push_back(save);
        }

        if (occ0 >= 0) {  // replace some env in rules
            std::string            dst_str = rule;
            std::string            sub_str = "[occ]";
            std::string            new_str = utils::concat("[", occ0, "]");
            std::string::size_type pos     = 0;
            while ((pos = dst_str.find(sub_str)) != std::string::npos) {
                dst_str.replace(pos, sub_str.length(), new_str);
            }
            rule = dst_str;
        }

        std::shared_ptr<RuleEvaluator> record_rule(  //
            new RuleEvaluator(rule, _dataset, mode, save, nsamp_ptr[0]));
        _ruleset->registerRules(record_rule);

        it++;
    }
}

Status& Kernel_Recorder::initializeKernel_impl(Status& stat) {
    bool not_parsed = _ruleset->getRules().size() == 0;
    if (not_parsed) parse();
    return stat;
}

Status& Kernel_Recorder::executeKernel_impl(Status& stat) {
    for (auto& irule : _ruleset->getRules()) { irule->calculateResult(isamp_ptr[0]); }
    return stat;
}

Status& Kernel_Recorder::finalizeKernel_impl(Status& stat) {
    for (auto& irule : _ruleset->getRules()) irule->collectResult();
    return stat;
}

};  // namespace PROJECT_NS

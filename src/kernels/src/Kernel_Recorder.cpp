#include "kids/Kernel_Recorder.h"

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

void Kernel_Recorder::token(Param::JSON& j) {
    std::string rule, mode = "average", save = "res.dat";
    if (j.is_string()) {
        rule = j.get<std::string>();
        // mode = "every";
        // save = "traj_{TRAJID}.dat";
    } else if (j.is_object()) {
        if (j.count("rule") == 1) {
            rule = j["rule"].get<std::string>();
        } else {
            throw kids_error("parser rule error");
        }
        if (j.count("mode") == 1) mode = j["mode"].get<std::string>();
        if (j.count("save") == 1) save = j["save"].get<std::string>();
    } else if (j.is_array()) {  // compile with old format @deprecated!
        std::string v0   = j[0].get<std::string>();
        auto        ipos = v0.find_first_of("#");
        if (ipos != std::string::npos) rule = v0.substr(ipos + 1, v0.size());
        if (ipos != std::string::npos) v0 = v0.substr(0, ipos);
        std::string v1 = j[1].get<std::string>();
        ipos           = v1.find_first_of("#");
        if (ipos != std::string::npos) rule += v1.substr(ipos + 1, v1.size());
        if (ipos != std::string::npos) v1 = v1.substr(0, ipos);
        rule += utils::concat("A(", v0, ",", v1, ")");
    } else {
        throw std::runtime_error("recorder parse error");
    }

    if (std::find(opened_files.begin(), opened_files.end(), save) == opened_files.end()) {
        std::shared_ptr<RuleEvaluator> record_time_rule(  //
            new RuleEvaluator("time(t{flowcontrol}:R, pertimeunit{flowcontrol}:R)", _dataset, "copy", save,
                              nsamp_ptr[0]));
        _ruleset->registerRules(record_time_rule);
        opened_files.push_back(save);
    }

    if (occ0 >= 0) {
        std::string            dst_str = rule;
        std::string            sub_str = "[occ]";
        std::string            new_str = utils::concat("[", occ0, "]");
        std::string::size_type pos     = 0;
        while ((pos = dst_str.find(sub_str)) != std::string::npos) { dst_str.replace(pos, sub_str.length(), new_str); }
        rule = dst_str;
    }

    std::shared_ptr<RuleEvaluator> record_rule(  //
        new RuleEvaluator(rule, _dataset, mode, save, nsamp_ptr[0]));
    _ruleset->registerRules(record_rule);
}

Status& Kernel_Recorder::initializeKernel_impl(Status& stat) {
    bool  not_parsed = _ruleset->getRules().size() == 0;
    auto& json       = *(_param->pjson());
    if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
        for (auto& j : (json["result"])) token(j);
    }
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

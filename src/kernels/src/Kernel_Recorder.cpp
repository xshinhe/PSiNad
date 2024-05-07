#include "kids/Kernel_Recorder.h"

#include "kids/Einsum.h"
#include "kids/RecordedRule.h"
#include "kids/RecorderIO.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Recorder::getName() { return "Kernel_Recorder"; }

int Kernel_Recorder::getType() const { return utils::hash(FUNCTION_NAME); }

Kernel_Recorder::~Kernel_Recorder(){};

void Kernel_Recorder::setInputParam_impl(std::shared_ptr<Param>& PM) {
    dt        = PM->get_double("dt", LOC(), phys::time_d);
    t0        = PM->get_double("t0", LOC(), phys::time_d, 0.0f);
    time_unit = PM->get_double("time_unit", LOC(), phys::time_d, 1.0f);
    directory = PM->get_string("directory", LOC(), "default");
}

void Kernel_Recorder::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    istep_ptr                     = DS->def(DATA::iter::istep);
    sstep_ptr                     = DS->def(DATA::iter::sstep);
    isamp_ptr                     = DS->def(DATA::iter::isamp);
    nsamp_ptr                     = DS->def(DATA::iter::nsamp);
    at_samplingstep_initially_ptr = DS->def(DATA::iter::at_samplingstep_initially);
}

void Kernel_Recorder::token(Param::JSON& j) {
    std::string rule, mode, save;
    if (j.is_string()) {
        rule = j.get<std::string>();
        mode = "average";
        save = "res.dat";
    } else if (j.is_object()) {
        if (j.count("rule") == 1) {
            rule = j["rule"].get<std::string>();
        } else {
            throw kids_error("parser rule error");
        }
        if (j.count("mode") == 1) mode = j["mode"].get<std::string>();
        if (j.count("save") == 1) save = j["save"].get<std::string>();
    } else if (j.is_array()) {
        std::string v0   = j[0].get<std::string>();
        auto        ipos = v0.find_first_of("#");
        if (ipos != std::string::npos) rule = v0.substr(ipos + 1, v0.size());
        if (ipos != std::string::npos) v0 = v0.substr(0, ipos);
        std::string v1 = j[1].get<std::string>();
        ipos           = v1.find_first_of("#");
        if (ipos != std::string::npos) rule += v1.substr(ipos + 1, v1.size());
        if (ipos != std::string::npos) v1 = v1.substr(0, ipos);
        rule += utils::concat("A(", v0, ",", v1, ")");
        mode = "average";
        save = "res.dat";
    } else {
        throw std::runtime_error("recorder parse error");
    }

    if (std::find(opened_files.begin(), opened_files.end(), save) == opened_files.end()) {
        std::shared_ptr<RecordedRule> record_time_rule(  //
            new RecordedRule("t{iter}", _dataset, "average", save, directory, nsamp_ptr[0]));
        Rules.push_back(record_time_rule);
        RecorderIO::registerRulesInRecorderIO(record_time_rule);
        opened_files.push_back(save);
    }

    std::shared_ptr<RecordedRule> record_rule(  //
        new RecordedRule(rule, _dataset, mode, save, directory, nsamp_ptr[0]));
    Rules.push_back(record_rule);
    RecorderIO::registerRulesInRecorderIO(record_rule);
}

Status& Kernel_Recorder::initializeKernel_impl(Status& stat) {
    bool  not_parsed = Rules.size() == 0;
    auto& json       = *(_param->pjson());
    if (not_parsed && json.count("result") == 1 && json["result"].is_array()) {
        for (auto& j : (json["result"])) token(j);
    }
    return stat;
}

Status& Kernel_Recorder::executeKernel_impl(Status& stat) {
    if (at_samplingstep_initially_ptr[0]) {
        for (auto& irecord : Rules) { irecord->calculateAtTimeSlice(isamp_ptr[0]); }
    }
    return stat;
}

};  // namespace PROJECT_NS

#include "kids/RuleSet.h"

#include <map>
#include <string>

#include "kids/RuleEvaluator.h"
#include "kids/concat.h"
#include "kids/fmt.h"

namespace PROJECT_NS {

std::vector<std::shared_ptr<RuleSet>>& RuleSet::getRuleSets() {
    static std::vector<std::shared_ptr<RuleSet>> _GLOBAL;
    return _GLOBAL;
}

void RuleSet::registerRulesInRuleSet(std::shared_ptr<RuleEvaluator>& exprRule) {
    auto& expression_ios = getRuleSets();
    for (auto it = expression_ios.begin(); it != expression_ios.end(); ++it) {
        if (exprRule->save == (*it)->unique_name) {
            (*it)->rules.push_back(exprRule);
            (*it)->appendHeader(exprRule);
            return;
        }
    }
    expression_ios.push_back(std::shared_ptr<RuleSet>(new RuleSet(exprRule)));
}

void RuleSet::registerRules(std::shared_ptr<RuleEvaluator>& exprRule) {
    rules.push_back(exprRule);
    RuleSet::registerRulesInRuleSet(exprRule);
}

RuleSet::RuleSet(std::shared_ptr<RuleEvaluator>& expr_rule) {
    unique_name      = expr_rule->save;
    header           = "";
    header_fstat     = "";
    totalFrameNumber = expr_rule->totalFrameNumber;
    rules.push_back(expr_rule);
    appendHeader(expr_rule);
}

void RuleSet::flush_all(const std::string& path, const std::string& suff, int level) {
    auto& expression_ios = getRuleSets();
    for (auto& io : expression_ios) io->flush(path, suff, level);
}

void RuleSet::flush(const std::string& path, const std::string& suff, int level) {
    std::ofstream ofs;
    if(level >=0) {
        ofs.open(utils::concat(path, "/", unique_name, suff));
        ofs << header << "\n";
    }else if(level==-2){
        ofs.open(utils::concat(path, "/fstat-", unique_name, suff));
        ofs << header_fstat << "\n";
    }else if(level==-1){
        ofs.open(utils::concat(path, "/", unique_name, suff), std::ios::app);
    }else{
        throw std::runtime_error("bad level");
    }
    for (int iframe = 0; iframe < totalFrameNumber; ++iframe) {
        for (auto& r : rules) {
            switch (level) {
                case -2: 
                case -1: 
                {
                    r->writeTo(ofs, r->result->dataPointerRes0, totalFrameNumber-1);
                    break;
                }
                case 0: {
                    r->writeTo(ofs, r->result->dataPointerRes0, iframe);
                    break;
                }
                case 1: {
                    r->writeTo(ofs, r->result->dataPointerRes1, iframe);
                    break;
                }
                case 2: {
                    r->writeTo(ofs, r->result->dataPointerRes2, iframe);
                    break;
                }
            }
        }
        ofs << "\n";
    }
    ofs.close();
}

void RuleSet::appendHeader(std::shared_ptr<RuleEvaluator>& expr_rule) {
    auto& res = expr_rule->result;
    if (!res->isTabular) return;

    std::stringstream ss;
    switch (res->dataType) {
        case kids_real_type: {
            for (int i = 0; i < expr_rule->totalTermNumber; ++i) {  //
                if (expr_rule->totalTermNumber == 1) {
                    ss << FMT(8) << res->name;
                } else {
                    ss << FMT(8) << utils::concat(res->name, "(", i, ")");
                }
            }
            break;
        }
        case kids_complex_type: {
            for (int i = 0; i < expr_rule->totalTermNumber; ++i) {
                if (expr_rule->totalTermNumber == 1) {
                    ss << FMT(8) << utils::concat("R:", res->name)  //
                       << FMT(8) << utils::concat("I:", res->name);
                } else {
                    ss << FMT(8) << utils::concat("R:", res->name, "(", i, ")")  //
                       << FMT(8) << utils::concat("I:", res->name, "(", i, ")");
                }
            }
            break;
        }
    }

    if(expr_rule->mode == "fstat") header_fstat += ss.str();
    header += ss.str();
}

std::vector<std::shared_ptr<RuleEvaluator>>& RuleSet::getRules() { return rules; }

Result RuleSet::getResult() {
    Result res{};
    for (auto& r : rules) {
        if (!r->result->isTabular) continue;
        res._data.push_back(std::make_tuple(r->result->name,                  //
                                            r->result->dataPointerRes0,       //
                                            r->result->dataType,              //
                                            r->result->stackedshape->size(),  //
                                            totalFrameNumber                  //
                                            ));
    }
    return res;
}

Result RuleSet::getCollect() {
    Result res{};
    for (auto& r : rules) {
        if (!r->result->isTabular) continue;
        res._data.push_back(std::make_tuple(r->result->name,                  //
                                            r->result->dataPointerRes1,       //
                                            r->result->dataType,              //
                                            r->result->stackedshape->size(),  //
                                            totalFrameNumber                  //
                                            ));
    }
    return res;
}

Result RuleSet::getReduced() {
    Result res{};
    for (auto& r : rules) {
        if (!r->result->isTabular) continue;
        res._data.push_back(std::make_tuple(r->result->name,                  //
                                            r->result->dataPointerRes2,       //
                                            r->result->dataType,              //
                                            r->result->stackedshape->size(),  //
                                            totalFrameNumber                  //
                                            ));
    }
    return res;
}

};  // namespace PROJECT_NS
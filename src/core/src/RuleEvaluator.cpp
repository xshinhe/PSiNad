#include "psnd/RuleEvaluator.h"

#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "psnd/Einsum.h"
#include "psnd/Expression.h"
#include "psnd/debug_utils.h"
#include "psnd/vars_list.h"

namespace PROJECT_NS {

/** custumized einsum-like function operation
 * @tparam          T               data type
 * @param[in]       EH              EinsumHelper object
 * @param[in]       data_inputs     vector of pointers of data of input tensors
 * @param[inout]    data_output     pointer stored data of output tensor
 */
template <typename T, typename Tout>
void einsum_fun(EinsumHelper&                  EH,           //
                Expression<T>&                 FP,           //
                const std::vector<void*>&      data_inputs,  //
                const std::vector<psnd_dtype>& data_idtype,  //
                Tout*                          data_output   //
) {
    auto& einsum_dims      = EH.einsum_dims;
    auto& einsum_iposes    = EH.einsum_iposes;
    auto& ipos_inputs      = EH.ipos_inputs;
    auto& ipos_inputs_init = EH.ipos_inputs_init;
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
    memcpy(ipos_inputs.data(), ipos_inputs_init.data(), total_tensor * sizeof(std::size_t));
    bool reset_zero = true;
    for (std::size_t iloop = 0, ipos_output = 0; iloop < total_loop; ++iloop) {
        if (reset_zero) data_output[ipos_output] = Tout(0);

        for (int iten = 0; iten < total_tensor; ++iten) {  //
            switch (data_idtype[iten]) {
                case psnd_int_type:
                    holder[iten] = cast_at<T, psnd_int>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case psnd_real_type:
                    holder[iten] = cast_at<T, psnd_real>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case psnd_complex_type: {
                    holder[iten] = cast_at<T, psnd_complex>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                }
            }
        }
        data_output[ipos_output] += cast<Tout>(FP.evaluate(holder.data()));

        std::size_t i = imax;
        while (++einsum_iposes[i] == einsum_dims[i] && i > imin) { einsum_iposes[i--] = 0; }

        for (int iten = 0; iten < total_tensor; ++iten)  //
            ipos_inputs[iten] += dh_inputs[iten].mapldims[i];
        ipos_output += dh_output_mapldims[i];
    }
}

RuleEvaluator::RuleEvaluator(const std::string&        rule,
                             std::shared_ptr<DataSet>& DS,  //
                             const std::string&        mode_str,
                             const std::string&        save,  //
                             std::size_t               totalFrameNumber)

    : rule{rule},
      mode{RuleEvaluatorPolicy::_from(mode_str)},
      save{save},
      totalFrameNumber{totalFrameNumber},
      numCollects{0} {
    _dataset = DS;  // for debug
    // here "=" is the separator
    auto        ipos     = rule.find("=");
    std::string vars_str = (ipos == std::string::npos) ? rule : rule.substr(0, ipos);
    expressionString     = (ipos == std::string::npos) ? "" : rule.substr(ipos + 1, rule.size());
    expressionString.erase(0, expressionString.find_first_not_of(" "));
    expressionString.erase(expressionString.find_last_not_of(" ") + 1);

    // here "(,)" are the separators
    ipos                 = vars_str.find("(");
    std::string res0_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(0, ipos);
    res0_str.erase(0, res0_str.find_first_not_of(" "));
    res0_str.erase(res0_str.find_last_not_of(" ") + 1);
    vars_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(ipos + 1, vars_str.size());
    char ch  = 'a';

    has_parameter        = vars_str.find("|") != std::string::npos;
    std::string vars_end = has_parameter ? "|" : ")";
    while ((ipos = vars_str.find(",")) != std::string::npos ||  //
           (ipos = vars_str.find(vars_end)) != std::string::npos) {
        std::string item_str = vars_str.substr(0, ipos);
        bool        get_end  = item_str.find("|") != std::string::npos;
        if (get_end) {
            ipos     = vars_str.find(vars_end);
            item_str = vars_str.substr(0, ipos);
        }
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str      = vars_str.substr(ipos + 1, vars_str.size());
        auto instance = VariableDescriptor(item_str, save, VariableDescriptorPolicy::Input);
        variables.push_back(std::move(instance));
        if (get_end) break;
    }
    if (variables.size() == 0) {
        std::string item_str = res0_str;
        auto        instance = VariableDescriptor(item_str, save, VariableDescriptorPolicy::Input);
        variables.push_back(std::move(instance));
    }
    for (auto& v : variables) v.defineIn(DS);

    if (has_parameter) {
        ipos = vars_str.find(",");
        if (ipos == std::string::npos) throw psnd_error(utils::concat("bad rule format : ", rule));
        std::string item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str = vars_str.substr(ipos + 1, vars_str.size());
        c1       = std::shared_ptr<VariableDescriptor>(
            new VariableDescriptor(item_str, save, VariableDescriptorPolicy::Input));

        ipos = vars_str.find(",");
        if (ipos == std::string::npos) throw psnd_error(utils::concat("bad rule format : ", rule));
        item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str = vars_str.substr(ipos + 1, vars_str.size());
        c2       = std::shared_ptr<VariableDescriptor>(
            new VariableDescriptor(item_str, save, VariableDescriptorPolicy::Input));

        ipos = vars_str.find(",");
        if (ipos != std::string::npos) throw psnd_error(utils::concat("bad rule format : ", rule));
        ipos     = vars_str.find(")");
        item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str = vars_str.substr(ipos + 1, vars_str.size());
        balance  = std::shared_ptr<VariableDescriptor>(
            new VariableDescriptor(item_str, save, VariableDescriptorPolicy::Input));

        c1->defineIn(DS);
        c2->defineIn(DS);
        balance->defineIn(DS);
    }

    if (mode != RuleEvaluatorPolicy::stat) {
        result = std::shared_ptr<VariableDescriptor>(
            new VariableDescriptor(res0_str, save, VariableDescriptorPolicy::TabularOutput));
    } else {
        result = std::shared_ptr<VariableDescriptor>(
            new VariableDescriptor(res0_str, save, VariableDescriptorPolicy::MetaOutput));
    }

    // built einsum information
    einsumString = "";
    for (int i = 0, shift = 0; i < variables.size(); ++i) {
        if (i != 0) einsumString += ",";
        if (variables[i].index == "") {  // update index if necessary
            for (int k = 0; k < variables[i].shape->rank(); ++k) { variables[i].index += (char) ((int) 'a' + shift++); }
        }
        einsumString += variables[i].index;
        inputShapes.push_back(variables[i].shape->dims());
        inputData.push_back(variables[i].dataPointerTrace);
        inputDataTypes.push_back(variables[i].dataType);
    }
    if (result->index != "") { einsumString += utils::concat("->", result->index); }
    einsumHelper = std::shared_ptr<EinsumHelper>(new EinsumHelper(einsumString, inputShapes));

    // define result in dataset
    expressionType = psnd_void_type;
    for (auto& v : variables) {
        if (v.dataType != psnd_int_type && v.dataType != psnd_real_type && v.dataType != psnd_complex_type)
            throw psnd_error("access bad type");
        if (v.dataType > expressionType) expressionType = v.dataType;
    }
    psnd_dtype res_type = (result->dataType == psnd_void_type) ? expressionType : result->dataType;

    totalTermNumber = einsumHelper->dh_output.total_size;

    // reset totalFramenumber
    if (mode == RuleEvaluatorPolicy::stat) totalFrameNumber = 1;
    result->defineIn(DS, res_type, einsumHelper->dh_output.dims, totalFrameNumber);

    // built expression
    if (expressionString == "") {                     // check if expression is trivial
        for (int i = 0; i < variables.size(); ++i) {  //
            expressionString += (i == 0) ? utils::concat("_", i) : utils::concat(" * _", i);
        }
    }
    // anyway, we define a corresponding expression (even if it is trivial)
    std::vector<std::string> varslist;
    for (int i = 0; i < variables.size(); ++i) varslist.push_back(utils::concat("_", i));
    switch (expressionType) {
        case psnd_real_type: {
            expressionId = Expression<psnd_real>::registerExpression(expressionString, varslist);
            break;
        }
        case psnd_complex_type: {
            expressionId = Expression<psnd_complex>::registerExpression(expressionString, varslist);
            break;
        }
    }
}

void RuleEvaluator::calculateResult(int sampleIndex, bool update) {
    sampleIndex         = sampleIndex % totalFrameNumber;
    size_t initialIndex = sampleIndex * totalTermNumber;

    for (auto& var : variables) var.checkTrace(sampleIndex, update);  // update Trace

    switch (expressionType) {
        case psnd_real_type: {
            auto&& eval = Expression<psnd_real>::getExpressions()[expressionId];
            if (result->dataType == psnd_real_type) {
                psnd_real* tracedata = (psnd_real*) (result->dataPointerTrace);
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, tracedata);
                if (has_parameter) {
                    psnd_real* c1data = (psnd_real*) c1->dataPointerRaw;
                    psnd_real* c2data = (psnd_real*) c2->dataPointerRaw;
                    psnd_real* bldata = (psnd_real*) balance->dataPointerRaw;
                    for (int i = 0; i < totalTermNumber; ++i)
                        tracedata[i] = c1data[0] * tracedata[i] - c2data[0] * bldata[i];
                }
                if (result->vtype == VariableDescriptorPolicy::TabularOutput) {
                    psnd_real* resdata = (psnd_real*) (result->dataPointerRes0) + initialIndex;
                    for (int i = 0; i < totalTermNumber; ++i) resdata[i] = tracedata[i];
                }
            }
            if (result->dataType == psnd_complex_type) {
                psnd_complex* tracedata = (psnd_complex*) (result->dataPointerTrace);
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, tracedata);
                if (has_parameter) {
                    psnd_real* c1data = (psnd_real*) c1->dataPointerRaw;
                    psnd_real* c2data = (psnd_real*) c2->dataPointerRaw;
                    psnd_real* bldata = (psnd_real*) balance->dataPointerRaw;
                    for (int i = 0; i < totalTermNumber; ++i)
                        tracedata[i] = c1data[0] * tracedata[i] - c2data[0] * bldata[i];
                }
                if (result->vtype == VariableDescriptorPolicy::TabularOutput) {
                    psnd_complex* resdata = (psnd_complex*) (result->dataPointerRes0) + initialIndex;
                    for (int i = 0; i < totalTermNumber; ++i) resdata[i] = tracedata[i];
                }
            }
            break;
        }
        case psnd_complex_type: {
            auto&& eval = Expression<psnd_complex>::getExpressions()[expressionId];
            if (result->dataType == psnd_real_type) {
                psnd_real* tracedata = (psnd_real*) (result->dataPointerTrace);
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, tracedata);
                if (has_parameter) {
                    psnd_real* c1data = (psnd_real*) c1->dataPointerRaw;
                    psnd_real* c2data = (psnd_real*) c2->dataPointerRaw;
                    psnd_real* bldata = (psnd_real*) balance->dataPointerRaw;
                    for (int i = 0; i < totalTermNumber; ++i)
                        tracedata[i] = c1data[0] * tracedata[i] - c2data[0] * bldata[i];
                }
                if (result->vtype == VariableDescriptorPolicy::TabularOutput) {
                    psnd_real* resdata = (psnd_real*) (result->dataPointerRes0) + initialIndex;
                    for (int i = 0; i < totalTermNumber; ++i) resdata[i] = tracedata[i];
                }
            }
            if (result->dataType == psnd_complex_type) {
                psnd_complex* tracedata = (psnd_complex*) (result->dataPointerTrace);
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, tracedata);
                if (has_parameter) {
                    psnd_real* c1data = (psnd_real*) c1->dataPointerRaw;
                    psnd_real* c2data = (psnd_real*) c2->dataPointerRaw;
                    psnd_real* bldata = (psnd_real*) balance->dataPointerRaw;
                    for (int i = 0; i < totalTermNumber; ++i)
                        tracedata[i] = c1data[0] * tracedata[i] - c2data[0] * bldata[i];
                }
                if (result->vtype == VariableDescriptorPolicy::TabularOutput) {
                    psnd_complex* resdata = (psnd_complex*) (result->dataPointerRes0) + initialIndex;
                    for (int i = 0; i < totalTermNumber; ++i) resdata[i] = tracedata[i];
                }
            }
            break;
        }
    }
}

void RuleEvaluator::collectResult() {
    if (result->vtype != VariableDescriptorPolicy::TabularOutput) return;
    switch (result->dataType) {
        case psnd_real_type: {
            psnd_real* fromdata = (psnd_real*) result->dataPointerRes0;
            psnd_real* todata   = (psnd_real*) result->dataPointerRes1;
            if (mode == RuleEvaluatorPolicy::copy || numCollects == 0) {
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] = fromdata[i];
            } else if (mode == RuleEvaluatorPolicy::sum) {
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] += fromdata[i];
            } else if (mode == RuleEvaluatorPolicy::average) {
                psnd_real k1 = numCollects / (psnd_real) (numCollects + 1);
                psnd_real k2 = 1.0e0 - k1;
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] = k1 * todata[i] + k2 * fromdata[i];
            }
            break;
        }
        case psnd_complex_type: {
            psnd_complex* fromdata = (psnd_complex*) result->dataPointerRes0;
            psnd_complex* todata   = (psnd_complex*) result->dataPointerRes1;
            if (mode == RuleEvaluatorPolicy::copy || numCollects == 0) {
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] = fromdata[i];
            } else if (mode == RuleEvaluatorPolicy::sum) {
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] += fromdata[i];
            } else if (mode == RuleEvaluatorPolicy::average) {
                psnd_real k1 = numCollects / (psnd_real) (numCollects + 1);
                psnd_real k2 = 1.0e0 - k1;
                for (int i = 0; i < result->stackedshape->size(); ++i) todata[i] = k1 * todata[i] + k2 * fromdata[i];
            }
            break;
        }
        default: {
            throw psnd_error("error type");
        }
    }
    numCollects++;
}

void RuleEvaluator::writeTo(std::ofstream& ofs, void* data, int sampleIndex) {
    if (result->vtype == VariableDescriptorPolicy::Input) return;
    switch (result->dataType) {
        case psnd_real_type: {
            psnd_real* resdata = (psnd_real*) (data) + sampleIndex * totalTermNumber;
            for (int i = 0; i < totalTermNumber; ++i) ofs << FMT(8) << resdata[i];
            break;
        }
        case psnd_complex_type: {
            psnd_complex* resdata = (psnd_complex*) (data) + sampleIndex * totalTermNumber;
            for (int i = 0; i < totalTermNumber; ++i) ofs << FMT(8) << real(resdata[i]) << FMT(8) << imag(resdata[i]);
            break;
        }
    }
}

};  // namespace PROJECT_NS
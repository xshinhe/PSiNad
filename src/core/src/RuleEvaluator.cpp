#include "kids/RuleEvaluator.h"

#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "kids/Einsum.h"
#include "kids/Expression.h"
#include "kids/debug_utils.h"

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
    bool reset_zero = true;
    for (std::size_t iloop = 0, ipos_output = 0; iloop < total_loop; ++iloop) {
        if (reset_zero) data_output[ipos_output] = Tout(0);

        for (int iten = 0; iten < total_tensor; ++iten) {  //
            switch (data_idtype[iten]) {
                case kids_int_type:
                    holder[iten] = cast_at<T, kids_int>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case kids_real_type:
                    holder[iten] = cast_at<T, kids_real>(data_inputs[iten], ipos_inputs[iten]);
                    break;
                case kids_complex_type: {
                    holder[iten] = cast_at<T, kids_complex>(data_inputs[iten], ipos_inputs[iten]);
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

RuleEvaluator::RuleEvaluator(const std::string& rule, std::shared_ptr<DataSet>& DS,  //
                             const std::string& mode, const std::string& save,       //
                             std::size_t totalFrameNumber)

    : rule{rule}, mode{mode}, save{save}, totalFrameNumber{totalFrameNumber}, numCollects{0} {
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
    result         = std::shared_ptr<VariableDescriptor>(new VariableDescriptor(res0_str));
    result->field  = "result";
    collect        = std::shared_ptr<VariableDescriptor>(new VariableDescriptor(res0_str));
    collect->field = "collect";
    reduced        = std::shared_ptr<VariableDescriptor>(new VariableDescriptor(res0_str));
    reduced->field = "reduced";

    vars_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(ipos + 1, vars_str.size());
    char ch  = 'a';
    while ((ipos = vars_str.find(",")) != std::string::npos ||  //
           (ipos = vars_str.find(")")) != std::string::npos) {
        std::string item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str      = vars_str.substr(ipos + 1, vars_str.size());
        auto instance = VariableDescriptor(item_str);
        instance.referIn(DS);
        variables.push_back(std::move(instance));
    }
    if (variables.size() == 0) {
        auto instance = VariableDescriptor(res0_str);
        instance.referIn(DS);
        variables.push_back(std::move(instance));
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
        inputData.push_back(variables[i].dataPointer);
        inputDataTypes.push_back(variables[i].dataType);
    }
    if (result->index != "") einsumString += utils::concat("->", result->index);
    einsumHelper = std::shared_ptr<EinsumHelper>(new EinsumHelper(einsumString, inputShapes));

    // define result in dataset
    expressionType = kids_void_type;
    for (auto& v : variables) {
        if (v.dataType != kids_int_type && v.dataType != kids_real_type && v.dataType != kids_complex_type)
            throw kids_error("access bad type");
        if (v.dataType > expressionType) expressionType = v.dataType;
    }
    kids_dtype res_type = (result->dataType == kids_void_type) ? expressionType : result->dataType;

    std::vector<std::size_t> res_cxxshape;
    res_cxxshape.push_back(totalFrameNumber);
    totalTermNumber = 1;
    for (auto& i : einsumHelper->dh_output.dims) {
        totalTermNumber *= i;
        res_cxxshape.push_back(i);
    }

    result->defineIn(DS, res_type, res_cxxshape);
    collect->defineIn(DS, res_type, res_cxxshape);
    reduced->defineIn(DS, res_type, res_cxxshape);

    // built expression
    if (expressionString == "") {  // check if expression is trivial
        expressionCheck = (variables.size() == 1) ? -1 : -2;
        for (int i = 0; i < variables.size(); ++i) {  //
            expressionString += (i == 0) ? utils::concat("_", i) : utils::concat(" * _", i);
        }
    } else {
        expressionCheck = 0;
    }
    // anyway, we define a corresponding expression (even if it is trivial)
    std::vector<std::string> varslist;
    for (int i = 0; i < variables.size(); ++i) varslist.push_back(utils::concat("_", i));
    switch (expressionType) {
        case kids_real_type: {
            expressionId = Expression<kids_real>::registerExpression(expressionString, varslist);
            break;
        }
        case kids_complex_type: {
            expressionId = Expression<kids_complex>::registerExpression(expressionString, varslist);
            break;
        }
    }

#ifdef LOCAL_DEBUG
    std::cout                                                 //
        << LOC() << "RuleEvaluator Data:\n"                   //
        << ".totalTermNumber = " << totalTermNumber << "\n"   /**< Number of terms. */
        << ".totalFrameNumber = " << totalFrameNumber << "\n" /**< Number of samples. */
        << ".rule = " << rule << "\n"                         /**< The expression rule. */
        << ".mode = " << mode << "\n"                         /**< The mode of evaluation. */
        << ".path = " << path << "\n"                         /**< Path to save results. */
        << ".save = " << save << "\n"                         /**< File name to save results. */
        << ".result = " << result
        << "\n" /**< Result of the expression. */
        // << ".variables = " << variables << "\n"            /**< Variables in the expression. */
        // << ".inputShapes = " << inputShapes << "\n"        /**< Shapes of input data. */
        // << ".inputData = " << inputData << "\n"            /**< Input data. */
        // << ".inputDataTypes = " << inputDataTypes << "\n"  /**< Data types of input data. */
        << ".expressionString" << expressionString
        << "\n" /**< String representation of the expression. */
        // << ".einsumHelper = " << einsumHelper << "\n"      /**< Shared pointer to EinsumHelper. */
        << ".einsumString = " << einsumString << "\n"     /**< String representation of expression type. */
        << ".expressionType = " << expressionType << "\n" /**< Type of the expression. */
        << ".expressionId = " << expressionId << "\n";    /**< ID of the expression. */

    for (int i = 0; i < inputDataTypes.size(); ++i) std::cout << inputDataTypes[i] << ",";
    std::cout << result->shape->size() << "\n";

#endif  // LOCAL_DEBUG
}

void RuleEvaluator::calculateResult(int sampleIndex) {
    sampleIndex         = sampleIndex % totalFrameNumber;
    size_t initialIndex = sampleIndex * totalTermNumber;

    // fast calculation for the special case: only copy data
    if (expressionCheck == -1) {
        if (result->dataType == kids_real_type) {
            kids_real* resdata = (kids_real*) (result->dataPointer) + initialIndex;
            kids_real* vardata = (kids_real*) (variables[0].dataPointer);
            for (int i = 0; i < totalTermNumber; ++i) resdata[i] = vardata[i];
        }
        if (result->dataType == kids_complex_type) {
            kids_complex* resdata = (kids_complex*) (result->dataPointer) + initialIndex;
            kids_complex* vardata = (kids_complex*) (variables[0].dataPointer);
            for (int i = 0; i < totalTermNumber; ++i) resdata[i] = vardata[i];
        }
        return;
    }

    // fast calculation for the special case: only use usual einsum
    if (expressionCheck == -2 && false) {  // not safe of reinterpret_cast!
        if (result->dataType == kids_real_type) {
            kids_real*              resdata = (kids_real*) (result->dataPointer) + initialIndex;
            std::vector<kids_real*> inputDataR;
            for (int i = 0; i < inputData.size(); ++i) inputDataR.push_back((kids_real*) inputData[i]);
            einsum<kids_real>(*einsumHelper, inputDataR, resdata);
        }
        if (result->dataType == kids_complex_type) {
            kids_complex*              resdata = (kids_complex*) (result->dataPointer) + initialIndex;
            std::vector<kids_complex*> inputDataC;
            for (int i = 0; i < inputData.size(); ++i) inputDataC.push_back((kids_complex*) inputData[i]);
            einsum<kids_complex>(*einsumHelper, inputDataC, resdata);
        }
        return;
    }

    switch (expressionType) {
        case kids_real_type: {
            auto&& eval = Expression<kids_real>::getExpressions()[expressionId];
            if (result->dataType == kids_real_type) {
                kids_real* resdata = (kids_real*) (result->dataPointer) + initialIndex;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            if (result->dataType == kids_complex_type) {
                kids_complex* resdata = (kids_complex*) (result->dataPointer) + initialIndex;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            break;
        }
        case kids_complex_type: {
            auto&& eval = Expression<kids_complex>::getExpressions()[expressionId];
            if (result->dataType == kids_real_type) {
                kids_real* resdata = (kids_real*) (result->dataPointer) + initialIndex;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            if (result->dataType == kids_complex_type) {
                kids_complex* resdata = (kids_complex*) (result->dataPointer) + initialIndex;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            break;
        }
    }
}

void RuleEvaluator::collectResult() {
    switch (result->dataType) {
        case kids_real_type: {
            kids_real* fromdata = (kids_real*) result->dataPointer;
            kids_real* todata   = (kids_real*) collect->dataPointer;
            if (mode == "copy" || numCollects == 0) {
                for (int i = 0; i < result->shape->size(); ++i) todata[i] = fromdata[i];
            } else if (mode == "sum") {
                for (int i = 0; i < result->shape->size(); ++i) todata[i] += fromdata[i];
            } else if (mode == "average") {
                kids_real k1 = numCollects / (kids_real)(numCollects + 1);
                kids_real k2 = 1.0e0 - k1;
                for (int i = 0; i < result->shape->size(); ++i) todata[i] = k1 * todata[i] + k2 * fromdata[i];
            }
            break;
        }
        case kids_complex_type: {
            kids_complex* fromdata = (kids_complex*) result->dataPointer;
            kids_complex* todata   = (kids_complex*) collect->dataPointer;
            if (mode == "copy" || numCollects == 0) {
                for (int i = 0; i < result->shape->size(); ++i) todata[i] = fromdata[i];
            } else if (mode == "sum") {
                for (int i = 0; i < result->shape->size(); ++i) todata[i] += fromdata[i];
            } else if (mode == "average") {
                kids_real k1 = numCollects / (kids_real)(numCollects + 1);
                kids_real k2 = 1.0e0 - k1;
                for (int i = 0; i < result->shape->size(); ++i) todata[i] = k1 * todata[i] + k2 * fromdata[i];
            }
            break;
        }
    }
    numCollects++;
}

void RuleEvaluator::writeTo(std::ofstream& ofs, void* data, int sampleIndex) {
    switch (result->dataType) {
        case kids_real_type: {
            kids_real* resdata = (kids_real*) (data) + sampleIndex * totalTermNumber;
            for (int i = 0; i < totalTermNumber; ++i) ofs << FMT(8) << resdata[i];
            break;
        }
        case kids_complex_type: {
            kids_complex* resdata = (kids_complex*) (data) + sampleIndex * totalTermNumber;
            for (int i = 0; i < totalTermNumber; ++i) ofs << FMT(8) << real(resdata[i]) << FMT(8) << imag(resdata[i]);
            break;
        }
    }
}

};  // namespace PROJECT_NS
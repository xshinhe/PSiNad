#include "kids/RecordedRule.h"

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

TokenVariable::TokenVariable(const std::string& token_string, std::shared_ptr<DataSet>& DS, char& eidxBegin,
                             bool initiallyRefered)
    : tokenString{token_string} {
    // pattern of "name{field@time}<index>:type";
    const std::regex pattern("([^{<:]*)(?:\\{([^@<]*)(?:@(.*?))?\\})?(?:<([^>]*)>)?(?::([^:]*))?");
    std::smatch      match;
    if (std::regex_match(tokenString, match, pattern)) {
        std::tie(name, field, time, index, type) = std::make_tuple(  //
            match[1], match[2], match[3], match[4], match[5]);
    } else {
        throw kids_error("cannot parser variable pattern");
    }
    if (field == "0") field = "init";
    if (field == "I") field = "integrator";
    if (field == "M") field = "model";
    if (field == "R") field = "result";
    if (field == "" && initiallyRefered) field = "integrator";
    if (field == "" && !initiallyRefered) field = "result";
    if (time == "" || time == "T") time = "t";
    if (initiallyRefered) this->referIn(DS, eidxBegin);
}

void TokenVariable::referIn(std::shared_ptr<DataSet>& DS, char& eidxBegin) {
    std::string key                        = utils::concat(field, ".", name);
    std::tie(dataType, dataPointer, shape) = DS->obtain(key);
    if (index == "") {
        for (int k = 0; k < shape->rank(); ++k, eidxBegin++) index += eidxBegin;
    }
#define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
    std::cout << LOC() << "TokenVariable Data:\n"
              << ".tokenString = " << tokenString << "\n"  //
              << ".name = " << name << "\n"                //
              << ".field = " << field << "\n"              //
              << ".index = " << index << "\n"              //
              << ".type = " << type << "\n"                //
              << ".time = " << time << "\n"                //
              << ".dataType = " << dataType << "\n"        //
              << ".dataPointer = " << dataPointer << "\n"  //
              << ".shape = " << ((shape) ? shape->to_string() : "") << "\n";
#endif  // LOCAL_DEBUG
}

RecordedRule::RecordedRule(const std::string& rule, std::shared_ptr<DataSet>& DS, const std::string& mode,
                           const std::string& save, const std::string& path, int numSamples)

    : rule{rule}, mode{mode}, save{save}, path{path}, numSamples{numSamples} {
    char        ch = 'a';
    std::string res0_str;
    std::string vars_str;
    std::string item_str;
    // here "=" is the separator
    auto ipos        = rule.find("=");
    vars_str         = (ipos == std::string::npos) ? rule : rule.substr(0, ipos);
    expressionString = (ipos == std::string::npos) ? "" : rule.substr(ipos + 1, rule.size());
    expressionString.erase(0, expressionString.find_first_not_of(" "));
    expressionString.erase(expressionString.find_last_not_of(" ") + 1);

    // here "(,)" are the separators
    ipos     = vars_str.find("(");
    res0_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(0, ipos);
    res0_str.erase(0, res0_str.find_first_not_of(" "));
    res0_str.erase(res0_str.find_last_not_of(" ") + 1);
    result        = std::shared_ptr<TokenVariable>(new TokenVariable(res0_str, DS, ch, false));
    result->field = "result";

    vars_str = (ipos == std::string::npos) ? vars_str : vars_str.substr(ipos + 1, vars_str.size());
    while ((ipos = vars_str.find(",")) != std::string::npos) {
        item_str = vars_str.substr(0, ipos);
        item_str.erase(0, item_str.find_first_not_of(" "));
        item_str.erase(item_str.find_last_not_of(" ") + 1);
        vars_str = vars_str.substr(ipos + 1, vars_str.size());
        variables.push_back(TokenVariable(item_str, DS, ch, true));
        // vars.push_back(VarItem(item_str, DS));
    }
    ipos     = vars_str.find(")");
    item_str = vars_str.substr(0, ipos);
    item_str.erase(0, item_str.find_first_not_of(" "));
    item_str.erase(item_str.find_last_not_of(" ") + 1);
    variables.push_back(TokenVariable(item_str, DS, ch, true));
    // vars.push_back(VarItem(item_str, DS));

    // pattern of "tokenRes (tokenVar1, tokenVar2, ...) = expr";
    // const std::regex pattern("([^=(]+)(?:\\(\\s*(?:,?\\s*([^,)]+))*\\s*\\))?(?:\\s*=\\s*([^=]+))?");
    // std::smatch      match;
    // if (std::regex_match(rule, match, pattern)) {
    //     for (int i = 0; i < match.size(); ++i) std::cout << match[i] << "\n";
    //     std::cout << match.size() << "\n";
    //     int numVariables = match.size() - 3;
    //     if (numVariables == 1 && match[2] == "") numVariables = 0;
    //     result        = std::shared_ptr<TokenVariable>(new TokenVariable(match[1], DS, ch, false));
    //     result->field = "result";
    //     for (int i = 0; i < numVariables; ++i) variables.push_back(TokenVariable(match[i + 2], DS, ch, true));
    //     if (numVariables == 0) { variables.push_back(TokenVariable(match[1], DS, ch, true)); }
    //     expressionString = match[numVariables + 2];
    // } else {
    //     throw kids_error("unknown pattern!");
    // }

    // built einsum information
    einsumString = "";
    for (int i = 0; i < variables.size(); ++i) {
        if (i != 0) einsumString += ",";
        einsumString += variables[i].index;
        inputShapes.push_back(variables[i].shape->dims());
        inputData.push_back(variables[i].dataPointer);
        inputDataTypes.push_back(variables[i].dataType);
    }
    if (result->index != "") einsumString += utils::concat("->", result->index);
    einsumHelper = std::shared_ptr<EinsumHelper>(new EinsumHelper(einsumString, inputShapes));

    // define result in dataset
    std::vector<std::size_t> shape_output_allslices;
    shape_output_allslices.push_back(numSamples);  // n samples of result
    numTerms = 1;
    for (auto&& i : einsumHelper->dh_output.dims) {
        numTerms *= i;
        shape_output_allslices.push_back(i);  // shape of result
    }

    expressionType = kids_void_type;
    for (auto&& v : variables) {
        if (v.dataType != kids_int_type && v.dataType != kids_real_type && v.dataType != kids_complex_type)
            throw std::runtime_error("access bad type");
        if (v.dataType > expressionType) expressionType = v.dataType;
    }
    switch (expressionType) {
        case kids_real_type:
            DS->def_real(utils::concat(result->field, ".", result->name),  //
                         shape_output_allslices, "CUSTOM DEFINED");
            break;
        case kids_complex_type:
            DS->def_complex(utils::concat(result->field, ".", result->name),  //
                            shape_output_allslices, "CUSTOM DEFINED");
            break;
    }
    result->referIn(DS, ch);  // defined and refered result data
    std::cout << LOC() << result->shape->size() << "\n";

    // built expression
    std::vector<std::string> varslist;
    for (int i = 0; i < variables.size(); ++i) varslist.push_back(utils::concat("_", i));
    if (expressionString == "") {
        for (int i = 0; i < variables.size(); ++i) {  //
            expressionString += (i == 0) ? utils::concat("_", i) : utils::concat(" * _", i);
        }
    }
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
    std::cout                                     //
        << LOC() << "RecordedRule Data:\n"        //
        << ".numTerms = " << numTerms << "\n"     /**< Number of terms. */
        << ".numSamples = " << numSamples << "\n" /**< Number of samples. */
        << ".rule = " << rule << "\n"             /**< The expression rule. */
        << ".mode = " << mode << "\n"             /**< The mode of evaluation. */
        << ".path = " << path << "\n"             /**< Path to save results. */
        << ".save = " << save << "\n"             /**< File name to save results. */
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

void RecordedRule::calculateAtTimeSlice(int sampleIndex) {
    sampleIndex = sampleIndex % numSamples;
    switch (expressionType) {
        case kids_real_type: {
            auto&& eval = Expression<kids_real>::getExpressions()[expressionId];
            if (result->dataType == kids_real_type) {
                kids_real* resdata = (kids_real*) (result->dataPointer) + sampleIndex * numTerms;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            if (result->dataType == kids_complex_type) {
                kids_complex* resdata = (kids_complex*) (result->dataPointer) + sampleIndex * numTerms;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            break;
        }
        case kids_complex_type: {
            auto&& eval = Expression<kids_complex>::getExpressions()[expressionId];
            if (result->dataType == kids_real_type) {
                kids_real* resdata = (kids_real*) (result->dataPointer) + sampleIndex * numTerms;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            if (result->dataType == kids_complex_type) {
                kids_complex* resdata = (kids_complex*) (result->dataPointer) + sampleIndex * numTerms;
                einsum_fun(*einsumHelper, eval, inputData, inputDataTypes, resdata);
            }
            break;
        }
    }
}

void RecordedRule::writeTo(std::ofstream& ofs, void* data, int sampleIndex) {
    switch (result->dataType) {
        case kids_real_type: {
            kids_real* resdata = (kids_real*) (data) + sampleIndex * numTerms;
            for (int i = 0; i < numTerms; ++i) ofs << FMT(8) << resdata[i];
            break;
        }
        case kids_complex_type: {
            kids_complex* resdata = (kids_complex*) (data) + sampleIndex * numTerms;
            for (int i = 0; i < numTerms; ++i) ofs << FMT(8) << real(resdata[i]) << FMT(8) << imag(resdata[i]);
            break;
        }
    }
}

};  // namespace PROJECT_NS
#include "kids/ExpressionWrapper.h"

#include "kids/Expression.h"

namespace PROJECT_NS {

// std::vector<ExpressionWrapper>& ExpressionWrapper::getExpressionWrappers() {
//     static std::vector<ExpressionWrapper> _GLOBAL;
//     return _GLOBAL;
// }

// int ExpressionWrapper::registerExpressionWrapper(const std::string& str, std::shared_ptr<DataSet>& DS,
//                                                  const std::string& field) {
//     auto& all_expressionwrappers = getExpressionWrappers();
//     int   idx                    = 0;
//     for (auto it = all_expressionwrappers.begin(); it != all_expressionwrappers.end(); ++it, ++idx) {
//         if (str == it->parsed_string && field == it->field) return idx;
//     }
//     all_expressionwrappers.push_back(ExpressionWrapper(str, DS, field));
//     return all_expressionwrappers.size() - 1;
// }

// int ExpressionWrapper::getSize() { return size; }

// kids_dtype ExpressionWrapper::getResultType() { return res_type; }

// std::string ExpressionWrapper::getName() { return unique_name; }

// ExpressionWrapper::ExpressionWrapper(const std::string& str, std::shared_ptr<DataSet>& DS, const std::string& field)
//     : field{field} {
//     auto ipos      = str.find("#");
//     bool use_alias = (ipos != std::string::npos);

//     std::string str_alias = "";
//     if (use_alias) {
//         str_alias = str.substr(ipos + 1, str.size());
//         str_alias.erase(0, str_alias.find_first_not_of(" "));
//         str_alias.erase(str_alias.find_last_not_of(" ") + 1);
//     }

//     parsed_string = str.substr(0, ipos);
//     parsed_string.erase(0, parsed_string.find_first_not_of(" "));
//     parsed_string.erase(parsed_string.find_last_not_of(" ") + 1);

//     ipos         = parsed_string.find("=");
//     expressionId = (ipos == std::string::npos) ? (-1) : 0;  // [single variable] or [function formula]

//     if (expressionId == -1) {  ///< [single variable]
//         auto inode = DS->node(utils::concat(field, ".", parsed_string));
//         res_type   = inode->type();
//         switch (res_type) {
//             case kids_int_type: {
//                 size = static_cast<Tensor<kids_int>*>(inode)->size();
//                 break;
//             }
//             case kids_real_type: {
//                 size = static_cast<Tensor<kids_real>*>(inode)->size();
//                 break;
//             }
//             case kids_complex_type: {
//                 size = static_cast<Tensor<kids_complex>*>(inode)->size();
//                 break;
//             }
//         }
//         variables_nodes.push_back(inode);
//         variables_types.push_back(res_type);
//         variables_sizes.push_back(size);
//         dims        = 1;
//         unique_name = (str_alias == "") ? parsed_string : str_alias;
//     } else {                                                              ///< [function formula]
//         parsed_declaration = parsed_string.substr(0, ipos);               // such as R(x,y)
//         parsed_expression  = parsed_string.substr(ipos + 1, str.size());  // such as x^2 + cos(y)

//         char type_char = parsed_declaration[0];  // from {'R', 'C'} for kids_real and kids_complex type

//         // remove (...) around the variables
//         parsed_varslist = parsed_declaration;
//         ipos            = parsed_varslist.find("(");
//         parsed_varslist = parsed_varslist.substr(ipos + 1, parsed_varslist.size());  // remove "("
//         ipos            = parsed_varslist.find(")");                                 //
//         parsed_varslist = parsed_varslist.substr(0, ipos);                           // remove ")"

//         // split variables into tokens by ',' delim
//         std::vector<std::string> tokens;
//         std::string              tokens_str = parsed_varslist;
//         ipos                                = tokens_str.find(',');
//         while (ipos != std::string::npos) {
//             tokens.push_back(tokens_str.substr(0, ipos));  // no-trim
//             tokens_str = tokens_str.substr(ipos + 1, tokens_str.size());
//             ipos       = tokens_str.find(',');
//         }
//         tokens.push_back(tokens_str);

//         // complete the function formula information
//         dims = tokens.size();
//         size = 1;
//         for (auto& var : tokens) {
//             int   id = regis_ExpressionWrapper(var, DS, field);
//             auto& f  = ExpressionWrapper::GLOBAL[id];
//             variables_nodes.push_back(f.variables_nodes[0]);
//             variables_types.push_back(f.variables_types[0]);
//             variables_sizes.push_back(f.variables_sizes[0]);
//             size *= f.variables_sizes[0];
//         }

//         switch (type_char) {
//             case 'R':
//                 res_type     = kids_real_type;
//                 expressionId = FPARSER<kids_real>::regis_FPARSER(parsed_expression, parsed_varslist);
//                 break;
//             case 'C':
//                 res_type     = kids_complex_type;
//                 expressionId = FPARSER<kids_complex>::regis_FPARSER(parsed_expression, parsed_varslist);
//                 break;
//         }
//         unique_name = (str_alias == "") ? utils::concat("[F", expressionId, "]") : str_alias;
//         std::cout << "Using Unique Name : " << unique_name << " = " << parsed_string << std::endl;
//     }

//     // calculate leading dimension
//     variables_ldims.resize(dims);
//     for (int i = dims - 1; i >= 0; --i) {
//         variables_ldims[i] = (i == dims - 1) ? 1 : variables_sizes[i + 1] * variables_ldims[i + 1];
//     }
// }

};  // namespace PROJECT_NS

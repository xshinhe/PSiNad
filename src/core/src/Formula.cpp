#include "kids/Formula.h"

namespace PROJECT_NS {

// template <typename T>
// int FPARSER<T>::regis_FPARSER(const std::string& str, const std::string& vars_str) {
//     int idx = 0;
//     for (auto it = GLOBAL.begin(); it != GLOBAL.end(); ++it, ++idx) {
//         if (str == it->expression_string) return idx;
//     }
//     GLOBAL.push_back(FPARSER<T>(str, vars_str));
//     return GLOBAL.size() - 1;
// }

// template <typename T>
// FPARSER<T>::FPARSER(const std::string& str, const std::string& vars_str)
//     : expression_string{str}, variables_string{vars_str} {
//     int res = fparser.Parse(expression_string, variables_string);
//     if (res != -1) throw std::runtime_error(fparser.ErrorMsg());
// }

template <typename T>
std::vector<FPARSER<T>> FPARSER<T>::GLOBAL;

////////////////////////////////////////////////////////////////////////////////////////////////

// int Formula::regis_Formula(const std::string& str, DataSet* DS, const std::string& field) {
//     int idx = 0;
//     for (auto it = GLOBAL.begin(); it != GLOBAL.end(); ++it, ++idx) {
//         if (str == it->parsed_string && field == it->field) return idx;
//     }
//     GLOBAL.push_back(Formula(str, DS, field));
//     return GLOBAL.size() - 1;
// }

// template <typename T>
// T Formula::eval(int idx) {
//     // [single variable]
//     if (FPARSER_ID == -1) return *((T*) variables_nodes[0]->data() + idx);

//     // [function formula]
//     std::vector<T> value_list(dims);
//     for (int i = 0, ID0 = 0; i < dims; ++i) {
//         ID0           = idx / variables_ldims[i];
//         idx           = idx % variables_ldims[i];
//         value_list[i] = *((T*) variables_nodes[i]->data() + ID0);
//     }
//     return FPARSER<T>::GLOBAL[FPARSER_ID].eval(value_list.data());
// }

// Formula::Formula(const std::string& str, DataSet* DS, const std::string& field) : parsed_string{str}, field{field} {
//     auto ipos  = str.find("=");
//     FPARSER_ID = (ipos == std::string::npos) ? (-1) : 0;  // [single variable] or [function formula]

//     if (FPARSER_ID == -1) {  ///< [single variable]
//         auto inode = std::get<3>(DS->info(utils::concat(field, ".", str)));
//         if (inode == nullptr) throw state_undefined_key_error(str);

//         res_type = inode->type();

//         variables_nodes.push_back(inode);
//         variables_types.push_back(inode->type());
//         variables_sizes.push_back(inode->size());

//         dims        = 1;
//         size        = inode->size();
//         unique_name = str;
//     } else {                                                    ///< [function formula]
//         parsed_declaration = str.substr(0, ipos);               // such as R(x,y)
//         parsed_expression  = str.substr(ipos + 1, str.size());  // such as x^2 + cos(y)

//         char type_char = parsed_declaration[0];  // from {'R', 'C'} for kids_real and kids_complex type

//         // remove (...) around the variables
//         parsed_varslist = parsed_declaration;
//         ipos            = parsed_varslist.find("(");
//         parsed_varslist = parsed_varslist.substr(ipos + 1, parsed_varslist.size());  // remove "("
//         ipos            = parsed_varslist.find(")");                                 //
//         parsed_varslist = parsed_varslist.substr(0, ipos);                           // remove ")"

//         // split variables into tokens by ',' delim
//         std::vector<std::string> tokens;
//         std::string tokens_str = parsed_varslist;
//         ipos                   = tokens_str.find(',');
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
//             int id  = regis_Formula(var, DS, field);
//             auto& f = Formula::GLOBAL[id];
//             variables_nodes.push_back(f.variables_nodes[0]);
//             variables_types.push_back(f.variables_types[0]);
//             variables_sizes.push_back(f.variables_sizes[0]);
//             size *= f.variables_sizes[0];
//         }

//         switch (type_char) {
//             case 'R':
//                 res_type   = DataSet::Type::Real;
//                 FPARSER_ID = FPARSER<kids_real>::regis_FPARSER(parsed_expression, parsed_varslist);
//                 break;
//             case 'C':
//                 res_type   = DataSet::Type::Complex;
//                 FPARSER_ID = FPARSER<kids_complex>::regis_FPARSER(parsed_expression, parsed_varslist);
//                 break;
//         }

//         unique_name = utils::concat("[F", FPARSER_ID, "]");
//         std::cout << "Using Unique Name : " << unique_name << " = " << parsed_string << std::endl;
//     }

//     // calculate leading dimension
//     variables_ldims.resize(dims);
//     for (int i = dims - 1; i >= 0; --i) {
//         variables_ldims[i] = (i == dims - 1) ? 1 : variables_sizes[i + 1] * variables_ldims[i + 1];
//     }
// }

std::vector<Formula> Formula::GLOBAL;

};  // namespace PROJECT_NS

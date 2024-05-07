#ifndef KIDS_ExpressionWrapper_H
#define KIDS_ExpressionWrapper_H

#include <string>
#include <vector>

#include "kids/DataSet.h"
#include "kids/concat.h"

namespace PROJECT_NS {

// class ExpressionWrapper {
//    public:
//     static std::vector<ExpressionWrapper>& getExpressionWrappers();

//     static int registerExpressionWrapper(const std::string& str, std::shared_ptr<DataSet>& DS,
//                                          const std::string& field);

//     // @deprecated
//     template <typename T>
//     T eval(int idx) {
//         // [single variable]
//         if (expressionId == -1) return *(static_cast<Tensor<T>*>(variables_nodes[0])->data() + idx);

//         // [function formula]
//         std::vector<T> value_list(dims);
//         for (int i = 0, ID0 = 0; i < dims; ++i) {
//             ID0 = idx / variables_ldims[i];
//             idx = idx % variables_ldims[i];
//             // value_list[i] = *((T*) variables_nodes[i]->data() + ID0);
//             switch (variables_types[i]) {
//                 case kids_int_type:
//                     value_list[i] =
//                         cast_at<T, kids_int>(static_cast<Tensor<kids_int>*>(variables_nodes[i])->data(), ID0);
//                     break;
//                 case kids_real_type:
//                     value_list[i] =
//                         cast_at<T, kids_real>(static_cast<Tensor<kids_real>*>(variables_nodes[i])->data(), ID0);
//                     break;
//                 case kids_complex_type:
//                     value_list[i] =
//                         cast_at<T, kids_complex>(static_cast<Tensor<kids_complex>*>(variables_nodes[i])->data(),
//                         ID0);
//                     break;
//             }
//         }
//         return FPARSER<T>::GLOBAL[expressionId].eval(value_list.data());
//     }

//     int getSize();

//     kids_dtype getResultType();

//     std::string getName();

//    private:
//     ExpressionWrapper(const std::string& str, std::shared_ptr<DataSet>& DS, const std::string& field);

//     std::string parsed_string;
//     std::string parsed_expression;
//     std::string parsed_declaration;
//     std::string parsed_varslist;
//     std::string field;
//     std::string unique_name;

//     kids_dtype              res_type = kids_void_type;
//     std::vector<Node*>      variables_nodes;
//     std::vector<kids_dtype> variables_types;
//     std::vector<size_t>     variables_sizes;
//     std::vector<size_t>     variables_ldims;

//     std::vector<void*> variables_datas;

//     size_t dims = 0;
//     size_t size = 0;
//     size_t expressionId;
// };

// std::vector<ExpressionWrapper> ExpressionWrapper::GLOBAL;


};  // namespace PROJECT_NS

#endif  //  KIDS_ExpressionWrapper_H
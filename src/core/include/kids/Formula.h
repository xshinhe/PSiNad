#ifndef FORMULA_H
#define FORMULA_H

#include <string>
#include <vector>

#include "fparser/fparser.hh"
#include "kids/DataSet.h"
#include "kids/concat.h"

namespace PROJECT_NS {

template <typename T>
class FPARSER {
   public:
    static std::vector<FPARSER<T>> GLOBAL;

    static int regis_FPARSER(const std::string& str, const std::string& vars_str) {
        int idx = 0;
        for (auto it = GLOBAL.begin(); it != GLOBAL.end(); ++it, ++idx) {
            if (str == it->expression_string) return idx;
        }
        GLOBAL.push_back(FPARSER<T>(str, vars_str));
        return GLOBAL.size() - 1;
    }

    inline T eval(T* data) { return fparser.Eval(data); }

   private:
    std::string expression_string;
    std::string variables_string;

    FunctionParserBase<T> fparser;

    FPARSER(const std::string& str, const std::string& vars_str) : expression_string{str}, variables_string{vars_str} {
        int res = fparser.Parse(expression_string, variables_string);
        if (res != -1) throw std::runtime_error(fparser.ErrorMsg());
    }
};
// template <typename T>
// std::vector<FPARSER<T>> FPARSER<T>::GLOBAL;

class Formula {
   public:
    static std::vector<Formula> GLOBAL;

    static int regis_Formula(const std::string& str, std::shared_ptr<DataSet>& DS, const std::string& field) {
        int idx = 0;
        for (auto it = GLOBAL.begin(); it != GLOBAL.end(); ++it, ++idx) {
            if (str == it->parsed_string && field == it->field) return idx;
        }
        GLOBAL.push_back(Formula(str, DS, field));
        return GLOBAL.size() - 1;
    }

    template <typename T>
    T eval(int idx) {
        // [single variable]
        if (FPARSER_ID == -1) return *(static_cast<Tensor<T>*>(variables_nodes[0])->data() + idx);

        // [function formula]
        std::vector<T> value_list(dims);
        for (int i = 0, ID0 = 0; i < dims; ++i) {
            ID0 = idx / variables_ldims[i];
            idx = idx % variables_ldims[i];
            // value_list[i] = *((T*) variables_nodes[i]->data() + ID0);
            switch (variables_types[i]) {
                case kids_int_type:
                    value_list[i] =
                        cast_at<T, kids_int>(static_cast<Tensor<kids_int>*>(variables_nodes[i])->data(), ID0);
                    break;
                case kids_real_type:
                    value_list[i] =
                        cast_at<T, kids_real>(static_cast<Tensor<kids_real>*>(variables_nodes[i])->data(), ID0);
                    break;
                case kids_complex_type:
                    value_list[i] =
                        cast_at<T, kids_complex>(static_cast<Tensor<kids_complex>*>(variables_nodes[i])->data(), ID0);
                    break;
            }
        }
        return FPARSER<T>::GLOBAL[FPARSER_ID].eval(value_list.data());
    }

    int get_size() { return size; }

    kids_dtype get_res_type() { return res_type; }

    std::string name() { return unique_name; }

   private:
    Formula(const std::string& str, std::shared_ptr<DataSet>& DS, const std::string& field);

    std::string parsed_string;
    std::string parsed_expression;
    std::string parsed_declaration;
    std::string parsed_varslist;
    std::string field;
    std::string unique_name;

    kids_dtype res_type = kids_void_type;
    std::vector<Node*> variables_nodes;
    std::vector<kids_dtype> variables_types;
    std::vector<int> variables_sizes;
    std::vector<int> variables_ldims;

    std::vector<void*> variables_datas;

    int dims = 0;
    int size = 0;
    int FPARSER_ID;
};

// std::vector<Formula> Formula::GLOBAL;


};  // namespace PROJECT_NS

#endif  //  FORMULA_H
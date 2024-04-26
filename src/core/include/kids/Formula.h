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
    Formula(const std::string& str, std::shared_ptr<DataSet>& DS, const std::string& field) : field{field} {
        auto ipos      = str.find("#");
        bool use_alias = (ipos != std::string::npos);

        std::string str_alias = "";
        if (use_alias) {
            str_alias = str.substr(ipos + 1, str.size());
            str_alias.erase(0, str_alias.find_first_not_of(" "));
            str_alias.erase(str_alias.find_last_not_of(" ") + 1);
        }

        parsed_string = str.substr(0, ipos);
        parsed_string.erase(0, parsed_string.find_first_not_of(" "));
        parsed_string.erase(parsed_string.find_last_not_of(" ") + 1);

        ipos       = parsed_string.find("=");
        FPARSER_ID = (ipos == std::string::npos) ? (-1) : 0;  // [single variable] or [function formula]

        if (FPARSER_ID == -1) {  ///< [single variable]
            auto inode = DS->node(utils::concat(field, ".", parsed_string));
            res_type   = inode->type();
            switch (res_type) {
                case kids_int_type: {
                    size = static_cast<Tensor<kids_int>*>(inode)->size();
                    break;
                }
                case kids_real_type: {
                    size = static_cast<Tensor<kids_real>*>(inode)->size();
                    break;
                }
                case kids_complex_type: {
                    size = static_cast<Tensor<kids_complex>*>(inode)->size();
                    break;
                }
            }
            variables_nodes.push_back(inode);
            variables_types.push_back(res_type);
            variables_sizes.push_back(size);
            dims        = 1;
            unique_name = (str_alias == "") ? parsed_string : str_alias;
        } else {                                                              ///< [function formula]
            parsed_declaration = parsed_string.substr(0, ipos);               // such as R(x,y)
            parsed_expression  = parsed_string.substr(ipos + 1, str.size());  // such as x^2 + cos(y)

            char type_char = parsed_declaration[0];  // from {'R', 'C'} for kids_real and kids_complex type

            // remove (...) around the variables
            parsed_varslist = parsed_declaration;
            ipos            = parsed_varslist.find("(");
            parsed_varslist = parsed_varslist.substr(ipos + 1, parsed_varslist.size());  // remove "("
            ipos            = parsed_varslist.find(")");                                 //
            parsed_varslist = parsed_varslist.substr(0, ipos);                           // remove ")"

            // split variables into tokens by ',' delim
            std::vector<std::string> tokens;
            std::string              tokens_str = parsed_varslist;
            ipos                                = tokens_str.find(',');
            while (ipos != std::string::npos) {
                tokens.push_back(tokens_str.substr(0, ipos));  // no-trim
                tokens_str = tokens_str.substr(ipos + 1, tokens_str.size());
                ipos       = tokens_str.find(',');
            }
            tokens.push_back(tokens_str);

            // complete the function formula information
            dims = tokens.size();
            size = 1;
            for (auto& var : tokens) {
                int   id = regis_Formula(var, DS, field);
                auto& f  = Formula::GLOBAL[id];
                variables_nodes.push_back(f.variables_nodes[0]);
                variables_types.push_back(f.variables_types[0]);
                variables_sizes.push_back(f.variables_sizes[0]);
                size *= f.variables_sizes[0];
            }

            switch (type_char) {
                case 'R':
                    res_type   = kids_real_type;
                    FPARSER_ID = FPARSER<kids_real>::regis_FPARSER(parsed_expression, parsed_varslist);
                    break;
                case 'C':
                    res_type   = kids_complex_type;
                    FPARSER_ID = FPARSER<kids_complex>::regis_FPARSER(parsed_expression, parsed_varslist);
                    break;
            }
            unique_name = (str_alias == "") ? utils::concat("[F", FPARSER_ID, "]") : str_alias;
            std::cout << "Using Unique Name : " << unique_name << " = " << parsed_string << std::endl;
        }

        // calculate leading dimension
        variables_ldims.resize(dims);
        for (int i = dims - 1; i >= 0; --i) {
            variables_ldims[i] = (i == dims - 1) ? 1 : variables_sizes[i + 1] * variables_ldims[i + 1];
        }
    }

    std::string parsed_string;
    std::string parsed_expression;
    std::string parsed_declaration;
    std::string parsed_varslist;
    std::string field;
    std::string unique_name;

    kids_dtype              res_type = kids_void_type;
    std::vector<Node*>      variables_nodes;
    std::vector<kids_dtype> variables_types;
    std::vector<int>        variables_sizes;
    std::vector<int>        variables_ldims;

    std::vector<void*> variables_datas;

    int dims = 0;
    int size = 0;
    int FPARSER_ID;
};

// std::vector<Formula> Formula::GLOBAL;


};  // namespace PROJECT_NS

#endif  //  FORMULA_H
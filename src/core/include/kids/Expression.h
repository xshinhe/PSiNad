#ifndef KIDS_EXPRESSION_H
#define KIDS_EXPRESSION_H

#include <string>
#include <vector>

#include "exprtk/exprtk.hpp"
#include "fparser/fparser.hh"
#include "kids/DataSet.h"
#include "kids/concat.h"

namespace PROJECT_NS {

template <typename T>
class Expression {
   public:
    static std::vector<Expression<T>>& getExpressions() {
        static std::vector<Expression<T>> _GLOBAL;
        return _GLOBAL;
    };

    static int registryExpression(const std::string& str, const std::string& vars_str) {
        auto& all_expressions = getExpressions();
        int idx               = 0;
        for (auto it = all_expressions.begin(); it != all_expressions.end(); ++it, ++idx) {
            if (str == it->expression_string) return idx;
        }
        all_expressions.push_back(Expression<T>(str, vars_str));
        return all_expressions.size() - 1;
    }

    // inline T eval(T* data) { return parser.Eval(data); }

    // private:
    std::string expression_string;
    std::string variables_string;
    // FunctionParserBase<T> parser;

    Expression(const std::string& str, const std::string& vars_str)
        : expression_string{str}, variables_string{vars_str} {
        // int res = parser.Parse(expression_string, variables_string);
        // if (res != -1) throw std::runtime_error(fparser.ErrorMsg());
    }
};

};  // namespace PROJECT_NS

// int main() {
//     // assert(1 > 2);
//     Expression<double>::registryExpression("x,y", "12,2");
//     Expression<double>::registryExpression("x,y,<d>", "12,2");
//     Expression<double>::registryExpression("x,y,<dd>", "12,2");
//     Expression<int>::registryExpression("x,y,x", "12,2");
//     Expression<int>::registryExpression("x,y,x,{i}", "12,2");
//     for (auto it = Expression<double>::getExpressions().begin(); it != Expression<double>::getExpressions().end();
//          ++it) {
//         std::cout << it->expression_string << "\n";
//     }
//     std::cout << "###\n";
//     for (auto it = Expression<int>::getExpressions().begin(); it != Expression<int>::getExpressions().end(); ++it) {
//         std::cout << it->expression_string << "\n";
//     }
//     return 0;
// }


#endif  //  KIDS_EXPRESSION_H

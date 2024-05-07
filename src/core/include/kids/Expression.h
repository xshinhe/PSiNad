/**@file        Expression.h
 * @brief       this file provide expression & parser
 * @details
 *
 * @author      Xin He
 * @date        2024-04
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-29  <td> moved from Formula.h
 * <tr><td> 2024-05-06  <td> improved & make clean
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_EXPRESSION_H
#define KIDS_EXPRESSION_H

#include <string>
#include <vector>

#ifdef KIDS_EXPRESSION_USE_EXPRTK
#include "exprtk/exprtk.hpp"               ///< @ref https://www.partow.net/programming/exprtk
#elif defined(KIDS_EXPRESSION_USE_LEPTON)  // LEPTON is a math parser by OPENMM-team
#else                                      // by default using fparser.h
#include "fparser/fparser.hh"              ///< @ref http://warp.povusers.org/FunctionParser/fparser.html
#endif                                     // KIDS_EXPRESSION_USE_EXPRTK

namespace PROJECT_NS {

/**
 * @brief A class representing a mathematical expression.
 *
 * This class allows registration, storage, and evaluation of mathematical expressions.
 * @tparam T The data type of the expression (e.g., float, double, complex)
 */
template <typename T>
class Expression {
   public:
    /**
     * @brief Get the list of registered expressions.
     *
     * @return A reference to the list of registered expressions.
     */
    static std::vector<Expression<T>>& getExpressions() {
        static std::vector<Expression<T>> _GLOBAL_EXPRESSIONS;
        return _GLOBAL_EXPRESSIONS;
    }

    /**
     * @brief Register a new mathematical expression.
     *
     * @param expression The mathematical expression string.
     * @param variables The list of variables in the expression.
     * @return The index of the registered expression.
     */
    static int registerExpression(const std::string& expression, const std::vector<std::string>& variables) {
        auto& allExpressions = getExpressions();
        int   index          = 0;
        for (auto it = allExpressions.begin(); it != allExpressions.end(); ++it, ++index) {
            if (expression == it->expressionString) return index;  // Check if expression already exists
        }
        allExpressions.push_back(Expression<T>(expression, variables));
        return allExpressions.size() - 1;
    }

    /**
     * @brief Evaluate the expression with given variable values.
     *
     * @param data An array containing values of variables.
     * @return The result of the expression evaluation.
     */
    inline T evaluate(T* data) { return _evaluate(data); }

   private:
#ifdef KIDS_EXPRESSION_USE_EXPRTK
    std::string              expressionString;  ///< The mathematical expression string.
    std::vector<std::string> variablesString;   ///< List of variables in the expression.
    exprtk::expression<T>    expressionValue;   ///< Expression object.
    std::vector<T>           variablesValue;    ///< Values of variables in the expression.

    /**
     * @brief Private constructor to initialize the expression.
     *
     * @param expressionStr The mathematical expression string.
     * @param varsStr The list of variables in the expression.
     */
    Expression(const std::string& expressionStr, const std::vector<std::string>& varsStr)
        : expressionString{expressionStr}, variablesString{varsStr} {
        exprtk::symbol_table<T> symbolTable;
        for (int i = 0; i < variablesString.size(); ++i)
            symbolTable.add_variable(variablesString[i], variablesValue[i]);
        symbolTable.add_constants();
        expressionValue.register_symbol_table(symbolTable);
        exprtk::parser<T>{}.compile(expressionString, expressionValue);
    }

    T _evaluate(T* data) {
        for (int i = 0; i < variablesValue.size(); ++i) variablesValue[i] = data[i];
        return expressionValue.value();
    }

#elif defined(KIDS_EXPRESSION_USE_LEPTON)  // LEPTON is a math parser by OPENMM-team
#else
    std::string           expressionString;  ///< The mathematical expression string.
    std::string           variablesString;   ///< List of variables in the expression.
    FunctionParserBase<T> parser;            ///< The parser for the expression.

    /**
     * @brief Private constructor to initialize the expression.
     *
     * @param expressionStr The mathematical expression string.
     * @param varsStr The list of variables in the expression.
     */
    Expression(const std::string& expressionStr, const std::vector<std::string>& varsStr)
        : expressionString{expressionStr} {
        variablesString = "";
        for (int i = 0; i < varsStr.size(); ++i) {
            variablesString += (i == 0) ? varsStr[i] : (std::string(",") + varsStr[i]);
        }
        int res = parser.Parse(expressionString, variablesString);
        if (res != -1) throw std::runtime_error(parser.ErrorMsg());
    }

    T _evaluate(T* data) { return parser.Eval(data); }

#endif  // KIDS_EXPRESSION_USE_EXPRTK
};

};  // namespace PROJECT_NS

#endif  //  KIDS_EXPRESSION_H

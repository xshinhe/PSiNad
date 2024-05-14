/**@file        VariableDescriptor.h
 * @brief       provide VariableDescriptor class
 * @details     VariableDescriptor
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
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-05-14  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_VariableDescriptor_H
#define KIDS_VariableDescriptor_H

#include "kids/DataSet.h"
#include "kids/Einsum.h"

namespace PROJECT_NS {

/**
 * @brief Represents a variable token in an expression rule.
 */
struct VariableDescriptor {
   public:
    /**
     * @brief Constructs a VariableDescriptor object.
     *
     * @param token_string The input string.
     * @param DS Shared pointer to the dataset.
     * @param initiallyDefined Flag indicating if the variable is initially defined (default is true).
     */
    VariableDescriptor(const std::string& token_string);

    /**
     * @brief Associates the variable with data in the DataSet.
     *
     * @param DS Shared pointer to the DataSet containing the variable's data.
     */
    void referIn(std::shared_ptr<DataSet>& DS);

    /**
     * @brief Define variable with data in the DataSet.
     *
     * @param DS Shared pointer to the DataSet containing the variable's data.
     */
    void defineIn(std::shared_ptr<DataSet>& DS, kids_dtype data_type, const std::vector<std::size_t>& cxxshape);

   private:
    friend class RuleEvaluator;
    friend class RuleSet;
    friend class Kernel;

    std::string tokenString; /**< The input token string. */
    std::string name;        /**< The name of the variable. */
    std::string field;       /**< The field of the variable. */
    std::string index;       /**< The index of the variable. */
    std::string type;        /**< The type of the variable. */
    std::string time;        /**< The time of the variable. */
    kids_dtype  dataType;    /**< The data type of the variable. */
    void*       dataPointer; /**< Pointer to the data of the variable. */
    Shape*      shape;       /**< Pointer to the shape of the variable. */
};

};  // namespace PROJECT_NS

#endif  // KIDS_VariableDescriptor_H

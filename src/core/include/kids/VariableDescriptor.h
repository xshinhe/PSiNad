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

struct RuleEvaluator;

/**
 * @brief Represents a variable token in an expression rule.
 */
struct VariableDescriptor {
   public:
    VariableDescriptor(const std::string& token_string, const std::string& save, bool is_output);

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
    void defineIn(std::shared_ptr<DataSet> DS, kids_dtype data_type = kids_void_type,
                  const std::vector<std::size_t>& cxxshape = {}, std::size_t totalFrameNumber = 0);

    void checkTrace(int sampleIndex);

   private:
    friend class RuleEvaluator;
    friend class RuleSet;
    friend class Kernel;

    bool isOutput;
    bool isTabular;

    // type: refer, def, tout

    std::string tokenString; /**< The input token string. */
    std::string field;       /**< The field of the variable. */
    std::string name;        /**< The name of the variable. */

    std::string save;  /**< The save of the variable. */
    std::string index; /**< The index of the variable. */
    std::string type;  /**< The type of the variable. */
    std::string time;  /**< The time of the variable. */

    std::string keyRaw;  // pointer to normal data
    std::string keyRec;  // transfer normal data in instant data during tracing
    std::string keyRes0;
    std::string keyRes1;
    std::string keyRes2;

    kids_dtype dataType;         /**< The data type of the variable. */
    void*      dataPointerRaw;   /**< Pointer to the data of the variable. */
    void*      dataPointerTrace; /**< Pointer to the data of the variable. */
    void*      dataPointerRes0;  /**< Pointer to the data of the variable. */
    void*      dataPointerRes1;  /**< Pointer to the data of the variable. */
    void*      dataPointerRes2;  /**< Pointer to the data of the variable. */
    Shape*     shape;            /**< Pointer to the shape of the variable. */
    Shape*     stackedshape;     /**< Pointer to the shape of the variable. */

    RuleEvaluator* RE;
};

};  // namespace PROJECT_NS

#endif  // KIDS_VariableDescriptor_H

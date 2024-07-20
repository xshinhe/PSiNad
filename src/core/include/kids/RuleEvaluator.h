/**@file        RuleEvaluator.h
 * @brief       provide RuleEvaluator class
 * @details     Build recording rule over DataSet
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

#ifndef KIDS_RuleEvaluator_H
#define KIDS_RuleEvaluator_H

#include "kids/DataSet.h"
#include "kids/Einsum.h"
#include "kids/VariableDescriptor.h"

namespace PROJECT_NS {

/**
 * @brief Represents a rule for evaluating an expression.
 */
struct RuleEvaluator {
    /**
     * @brief Constructs an RuleEvaluator object.
     *
     * @param rule The expression rule.
     * @param DS Shared pointer to the dataset.
     * @param mode The mode of evaluation (default is "average").
     * @param save The filename to save results (default is "res.dat").
     * @param nsamples Number of samples (default is 1).
     */
    RuleEvaluator(const std::string& rule, std::shared_ptr<DataSet>& DS,  //
                  const std::string& mode, const std::string& save,       //
                  std::size_t totalFrameNumber);

    /**
     * @brief Calculates the expression at a particular time slice.
     *
     * @param sampleIndex Index of the sample (i.e. frame index, default is 0).
     */
    void calculateResult(int sampleIndex = 0);

    void collectResult();

    /**
     * @brief Writes the results to a file stream.
     *
     * @param ofs The output file stream.
     * @param sampleIndex Index of the sample.
     */
    void writeTo(std::ofstream& ofs, void* data, int sampleIndex);

    std::size_t                           totalTermNumber;  /**< Number of terms. */
    std::size_t                           totalFrameNumber; /**< Number of samples. */
    std::size_t                           numCollects;      /**< Number of Collects. */
    std::string                           rule;             /**< The expression rule. */
    std::string                           mode;             /**< The mode of evaluation. */
    std::string                           save;             /**< File name to save results. */
    std::shared_ptr<VariableDescriptor>   result;           /**< Result of the expression. */
    std::vector<VariableDescriptor>       variables;        /**< Variables in the expression. */
    std::vector<std::vector<std::size_t>> inputShapes;      /**< Shapes of input data. */
    std::vector<void*>                    inputData;        /**< Input data. */
    std::vector<kids_dtype>               inputDataTypes;   /**< Data types of input data. */
    std::string                           expressionString; /**< String representation of the expression. */
    std::shared_ptr<EinsumHelper>         einsumHelper;     /**< Shared pointer to EinsumHelper. */
    std::string                           einsumString;     /**< String representation of expression type. */
    kids_dtype                            expressionType;   /**< Type of the expression. */
    size_t                                expressionId;     /**< ID of the expression. */

    bool                                has_parameter;
    std::shared_ptr<VariableDescriptor> c1;      /**< Result of the expression. */
    std::shared_ptr<VariableDescriptor> c2;      /**< Result of the expression. */
    std::shared_ptr<VariableDescriptor> balance; /**< Result of the expression. */
};

};  // namespace PROJECT_NS

#endif  // KIDS_RuleEvaluator_H

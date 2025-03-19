/**@file        RuleSet.h
 * @brief       provide RuleSet class
 * @details     A set of recording rules, and utils for wrapping io/stream
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

#ifndef KIDS_RULESET_H
#define KIDS_RULESET_H

#include <map>
#include <string>

#include "kids/Result.h"
#include "kids/RuleEvaluator.h"
#include "kids/concat.h"
#include "kids/fmt.h"

namespace PROJECT_NS {

/**
 * @brief Represents a handler for input/output operations related to expressions.
 *
 * This class manages the input/output operations for expressions. It allows registration of rules
 * and flushing the data to files. Each rule is associated with specific input/output file names.
 */
class RuleSet final {
   public:
    /**
     * @brief Get the list of RuleSet instances.
     *
     * @return A reference to the list of RuleSet instances.
     */
    static std::vector<std::shared_ptr<RuleSet>>& getRuleSets();

    /**
     * @brief Register rules in the RuleSet for a specific input/output file.
     *
     * If the RuleSet instance for the given file name already exists, the rule is added to it.
     * Otherwise, a new RuleSet instance is created.
     *
     * @param expr_rule Pointer to the RuleEvaluator to be registered.
     */
    static void registerRulesInRuleSet(std::shared_ptr<RuleEvaluator>& exprRule);

    /**
     * @brief Destructor.
     *
     * Closes the output file stream.
     */
    virtual ~RuleSet(){};

    /**
     * @brief Flush data to all output files.
     */
    static void flush_all(const std::string& path, const std::string& suff, int level);

    std::vector<std::shared_ptr<RuleEvaluator>>& getRules();

    Result getResult();

    Result getCollect();

    Result getReduced();

   private:
    friend class Kernel_Recorder;
    std::string                                 unique_name;
    std::string                                 header;            ///< The header of set
    std::string                                 header_fstat;
    size_t                                      totalFrameNumber;  ///< The number of frames.
    std::vector<std::shared_ptr<RuleEvaluator>> rules;             ///< List of RuleEvaluator
    bool                                        writeable = true;

    RuleSet() : writeable{false} {};  // @not registered

    void registerRules(std::shared_ptr<RuleEvaluator>& exprRule);

    /**
     * @brief Constructor.
     *
     * Initializes the RuleSet instance with the given RuleEvaluator and file name.
     * Opens the output file stream.
     *
     * @param expr_rule Pointer to the RuleEvaluator.
     */
    RuleSet(std::shared_ptr<RuleEvaluator>& expr_rule);

    /**
     * @brief Flush data to the output file.
     *
     * @param path Path for writing the file
     * Writes the file header followed by data for each frame according to registered rules.
     */
    void flush(const std::string& path, const std::string& suff, int level);

    /**
     * @brief Append header information to the output file.
     *
     * @param expr_rule Pointer to the RuleEvaluator.
     */
    void appendHeader(std::shared_ptr<RuleEvaluator>& expr_rule);
};
};  // namespace PROJECT_NS

#endif  // KIDS_RULESET_H

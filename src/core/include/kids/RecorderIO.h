#ifndef KIDS_EXPRESSIONIO_H
#define KIDS_EXPRESSIONIO_H

#include <map>
#include <string>

#include "kids/RecordedRule.h"
#include "kids/concat.h"
#include "kids/fmt.h"

namespace PROJECT_NS {

/**
 * @brief Represents a handler for input/output operations related to expressions.
 *
 * This class manages the input/output operations for expressions, allowing registration of rules and flushing the data
 * to files. Rules are associated with specific input/output file names.
 */
class RecorderIO final {
   public:
    /**
     * @brief Get RecorderIO instances.
     *
     * @return A reference to the dictionary of RecorderIO instances.
     */
    static std::vector<std::shared_ptr<RecorderIO>>& getRecorderIOs();

    /**
     * @brief Register rules in the RecorderIO for a specific input/output file.
     *
     * If the RecorderIO instance for the given file name already exists, the rule is added to it.
     * Otherwise, a new RecorderIO instance is created.
     *
     * @param expr_rule Pointer to the RecordedRule to be registered.
     * @param io_save The name of the input/output file.
     */
    static void registerRulesInRecorderIO(std::shared_ptr<RecordedRule>& exprRule) {
        auto& expression_ios = getRecorderIOs();
        for (auto it = expression_ios.begin(); it != expression_ios.end(); ++it) {
            if (exprRule->save == (*it)->save) {
                (*it)->rules.push_back(exprRule);
                (*it)->appendHeader(exprRule);
                return;
            }
        }
        expression_ios.push_back(std::shared_ptr<RecorderIO>(new RecorderIO(exprRule)));
    }

    /**
     * @brief Destructor.
     *
     * Closes the output file stream.
     */
    ~RecorderIO() { ofs.close(); }

    static void flush_all() {
        auto& expression_ios = getRecorderIOs();
        for (auto& io : expression_ios) io->flush();
    }

   private:
    /**
     * @brief Constructor.
     *
     * Initializes the RecorderIO instance with the given RecordedRule and file name.
     * Opens the output file stream.
     *
     * @param expr_rule Pointer to the RecordedRule.
     */
    RecorderIO(std::shared_ptr<RecordedRule>& expr_rule) : save{expr_rule->save}, header{""} {
        ofs.open(utils::concat(expr_rule->path, "/", expr_rule->save));
        nframe = expr_rule->numSamples;
        rules.push_back(expr_rule);
        appendHeader(expr_rule);
    }

    /**
     * @brief Flushes the data to the output file.
     *
     * Writes the file header followed by data for each frame according to registered rules.
     */
    void flush() {
        ofs << header << "\n";
        for (int iframe = 0; iframe < nframe; ++iframe) {
            for (auto& r : rules) { r->writeTo(ofs, r->result->dataPointer, iframe); }
            ofs << "\n";
        }
    }

    void appendHeader(std::shared_ptr<RecordedRule>& expr_rule) {
        auto&             res = expr_rule->result;
        std::stringstream ss;
        switch (res->dataType) {
            case kids_real_type: {
                for (int i = 0; i < expr_rule->numTerms; ++i) { ss << FMT(8) << utils::concat(res->name, i); }
                break;
            }
            case kids_complex_type: {
                for (int i = 0; i < expr_rule->numTerms; ++i) {
                    ss << FMT(8) << utils::concat("R", res->name, i)  //
                       << FMT(8) << utils::concat("I", res->name, i);
                }
                break;
            }
        }
        header += ss.str();
    }
    std::string                                save;
    std::string                                header;  ///< The header information for the file.
    size_t                                     nframe;  ///< The number of frames.
    std::ofstream                              ofs;     ///< Output file stream.
    std::vector<std::shared_ptr<RecordedRule>> rules;   ///< List of RecordedRule pointers
};

};  // namespace PROJECT_NS

#endif  // KIDS_EXPRESSIONIO_H

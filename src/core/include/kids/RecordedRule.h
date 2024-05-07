#ifndef KIDS_RecordedRule_H
#define KIDS_RecordedRule_H

#include "kids/DataSet.h"
#include "kids/Einsum.h"

namespace PROJECT_NS {

/**
 * @brief Represents a variable token in an expression rule.
 */
struct TokenVariable {
   public:
    /**
     * @brief Constructs a TokenVariable object.
     *
     * @param token_string The input string.
     * @param DS Shared pointer to the dataset.
     * @param initiallyDefined Flag indicating if the variable is initially defined (default is true).
     */
    TokenVariable(const std::string& token_string, std::shared_ptr<DataSet>& DS, char& eidxBegin,
                  bool initiallyRefered = true);

    /**
     * @brief Associates the variable with data in the DataSet.
     *
     * @param DS Shared pointer to the DataSet containing the variable's data.
     */
    void referIn(std::shared_ptr<DataSet>& DS, char& eidxBegin);

   private:
    friend class RecordedRule;
    friend class RecorderIO;
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

/**
 * @brief Represents a rule for evaluating an expression.
 */
struct RecordedRule {
    /**
     * @brief Constructs an RecordedRule object.
     *
     * @param rule The expression rule.
     * @param DS Shared pointer to the dataset.
     * @param mode The mode of evaluation (default is "average").
     * @param save The filename to save results (default is "res.dat").
     * @param path The path to save file (default is "default").
     * @param nsamples Number of samples (default is 1).
     */
    RecordedRule(const std::string& rule, std::shared_ptr<DataSet>& DS,  //
                 const std::string& mode = "average", const std::string& save = "res.dat",
                 const std::string& path = "default", int nsamples = 1);

    /**
     * @brief Calculates the expression at a particular time slice.
     *
     * @param sampleIndex Index of the sample (default is 0).
     */
    void calculateAtTimeSlice(int sampleIndex = 0);

    /**
     * @brief Writes the results to a file stream.
     *
     * @param ofs The output file stream.
     * @param sampleIndex Index of the sample.
     */
    void writeTo(std::ofstream& ofs, void* data, int sampleIndex);

    int                                   numTerms;         /**< Number of terms. */
    int                                   numSamples;       /**< Number of samples. */
    std::string                           rule;             /**< The expression rule. */
    std::string                           mode;             /**< The mode of evaluation. */
    std::string                           path;             /**< Path to save results. */
    std::string                           save;             /**< File name to save results. */
    std::shared_ptr<TokenVariable>        result;           /**< Result of the expression. */
    std::vector<TokenVariable>            variables;        /**< Variables in the expression. */
    std::vector<std::vector<std::size_t>> inputShapes;      /**< Shapes of input data. */
    std::vector<void*>                    inputData;        /**< Input data. */
    std::vector<kids_dtype>               inputDataTypes;   /**< Data types of input data. */
    std::string                           expressionString; /**< String representation of the expression. */
    std::shared_ptr<EinsumHelper>         einsumHelper;     /**< Shared pointer to EinsumHelper. */
    std::string                           einsumString;     /**< String representation of expression type. */
    kids_dtype                            expressionType;   /**< Type of the expression. */
    std::size_t                           expressionId;     /**< ID of the expression. */
};

};  // namespace PROJECT_NS

#endif  // KIDS_RecordedRule_H

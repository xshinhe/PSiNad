/** @file Policy.h
 *  @brief Provides macros for defining policy namespaces with enum-to-string mapping utilities
 *  @details This header defines a set of macros to simplify the creation of "policy" namespaces,
 *           each containing an enumeration, a string-to-enum mapping, and utility functions for
 *           conversion and documentation. Policies are used to encapsulate sets of valid options
 *           (e.g., algorithm types, parameter modes) with built-in validation and user assistance.
 *
 *  @author Xin He
 *  @date 2024-09
 */
#ifndef PSND_POLICY_H
#define PSND_POLICY_H

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

/**
 * @brief Converts a macro argument to a string literal
 * @param x Macro argument to convert
 * @return String literal representation of x
 */
#define VAR_NAME(x) #x

/**
 * @brief Concatenates two macro arguments into a single identifier
 * @param x First argument
 * @param y Second argument
 * @return Concatenated identifier x##y
 */
#define CAT_NAME(x, y) x##y

/**
 * @brief Selects a name by appending a number (e.g., NAME_5 for NUM=5)
 * @param NAME Base name
 * @param NUM Suffix number
 * @return NAME_NUM identifier
 */
#define SELECT_NAME(NAME, NUM) CAT_NAME(NAME##_, NUM)

/**
 * @brief Counts the number of variable arguments (up to 50)
 * @details Uses recursive macro expansion to determine argument count
 * @param ... Variable arguments to count
 * @return Number of arguments (0-50)
 */
#define ARG_COUNT(...)                                                  \
    ARG_COUNT_PRIVATE_IMPL(0, ##__VA_ARGS__,                       /**/ \
                           50, 49, 48, 47, 46, 45, 44, 43, 42, 41, /**/ \
                           40, 39, 38, 37, 36, 35, 34, 33, 32, 31, /**/ \
                           30, 29, 28, 27, 26, 25, 24, 23, 22, 21, /**/ \
                           20, 19, 18, 17, 16, 15, 14, 13, 12, 11, /**/ \
                           10, 9, 8, 7, 6, 5, 4, 3, 2, 1,          /**/ \
                           0)

/**
 * @brief Private implementation helper for ARG_COUNT
 * @details Uses variadic macro tricks to extract argument count
 * @param _0 Dummy parameter
 * @param ... Count values (50 down to 0)
 * @param count Selected count value (output)
 * @return count parameter
 */
#define ARG_COUNT_PRIVATE_IMPL(_0,                                               /**/ \
                               _1, _2, _3, _4, _5, _6, _7, _8, _9, _10,          /**/ \
                               _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, /**/ \
                               _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, /**/ \
                               _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, /**/ \
                               _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, /**/ \
                               count, ...)                                            \
    count

/**
 * @brief Selects the appropriate KV_TERMS_N macro based on argument count
 * @param NAME Base name ("KV_TERMS")
 * @param ... Variable arguments to process
 * @return Expanded KV_TERMS_N macro
 */
#define VA_SELECT(NAME, ...) SELECT_NAME(NAME, ARG_COUNT(__VA_ARGS__))(__VA_ARGS__)

/**
 * @brief Generates key-value pairs for policy enumerators (50 arguments)
 * @param TERM1 First enumerator
 * @param ... Remaining enumerators
 * @return Initializer list entries for map construction
 */
#define KV_TERMS_50(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_49(__VA_ARGS__)
#define KV_TERMS_49(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_48(__VA_ARGS__)
#define KV_TERMS_48(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_47(__VA_ARGS__)
#define KV_TERMS_47(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_46(__VA_ARGS__)
#define KV_TERMS_46(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_45(__VA_ARGS__)
#define KV_TERMS_45(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_44(__VA_ARGS__)
#define KV_TERMS_44(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_43(__VA_ARGS__)
#define KV_TERMS_43(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_42(__VA_ARGS__)
#define KV_TERMS_42(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_41(__VA_ARGS__)
#define KV_TERMS_41(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_40(__VA_ARGS__)
#define KV_TERMS_40(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_39(__VA_ARGS__)
#define KV_TERMS_39(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_38(__VA_ARGS__)
#define KV_TERMS_38(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_37(__VA_ARGS__)
#define KV_TERMS_37(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_36(__VA_ARGS__)
#define KV_TERMS_36(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_35(__VA_ARGS__)
#define KV_TERMS_35(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_34(__VA_ARGS__)
#define KV_TERMS_34(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_33(__VA_ARGS__)
#define KV_TERMS_33(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_32(__VA_ARGS__)
#define KV_TERMS_32(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_31(__VA_ARGS__)
#define KV_TERMS_31(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_30(__VA_ARGS__)
#define KV_TERMS_30(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_29(__VA_ARGS__)
#define KV_TERMS_29(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_28(__VA_ARGS__)
#define KV_TERMS_28(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_27(__VA_ARGS__)
#define KV_TERMS_27(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_26(__VA_ARGS__)
#define KV_TERMS_26(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_25(__VA_ARGS__)
#define KV_TERMS_25(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_24(__VA_ARGS__)
#define KV_TERMS_24(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_23(__VA_ARGS__)
#define KV_TERMS_23(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_22(__VA_ARGS__)
#define KV_TERMS_22(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_21(__VA_ARGS__)
#define KV_TERMS_21(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_20(__VA_ARGS__)
#define KV_TERMS_20(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_19(__VA_ARGS__)
#define KV_TERMS_19(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_18(__VA_ARGS__)
#define KV_TERMS_18(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_17(__VA_ARGS__)
#define KV_TERMS_17(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_16(__VA_ARGS__)
#define KV_TERMS_16(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_15(__VA_ARGS__)
#define KV_TERMS_15(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_14(__VA_ARGS__)
#define KV_TERMS_14(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_13(__VA_ARGS__)
#define KV_TERMS_13(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_12(__VA_ARGS__)
#define KV_TERMS_12(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_11(__VA_ARGS__)
#define KV_TERMS_11(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_10(__VA_ARGS__)
#define KV_TERMS_10(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_9(__VA_ARGS__)
#define KV_TERMS_9(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_8(__VA_ARGS__)
#define KV_TERMS_8(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_7(__VA_ARGS__)
#define KV_TERMS_7(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_6(__VA_ARGS__)
#define KV_TERMS_6(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_5(__VA_ARGS__)
#define KV_TERMS_5(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_4(__VA_ARGS__)
#define KV_TERMS_4(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_3(__VA_ARGS__)
#define KV_TERMS_3(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_2(__VA_ARGS__)
#define KV_TERMS_2(TERM1, ...) {VAR_NAME(TERM1), TERM1}, KV_TERMS_1(__VA_ARGS__)
#define KV_TERMS_1(TERM1, ...) {VAR_NAME(TERM1), TERM1}
#define KV_TERMS_0() /* Empty for 0 arguments */

/**
 * @brief Generates key-value pairs for policy enumerators
 * @details Expands to a list of {string, enum} pairs for map initialization
 * @param ... Policy enumerator values
 * @return Comma-separated initializer list entries
 */
#define KV_TERMS(...) VA_SELECT(KV_TERMS, ##__VA_ARGS__)


/**
 * @brief Defines a policy type with enumerators and helper functions
 * @details Generates:
 * - An enum class with specified enumerators + INVALID
 * - A map from string names to enum values
 * - A help function to print available values
 * - A conversion function from string to enum
 * @param Policy Name of the policy type
 * @param ... Enumerator values for the policy
 */
#define DEFINE_POLICY(Policy, ...)                                                    \
    namespace Policy {                                                                \
    /** @brief Enum type representing policy values */                                \
    enum _type { INVALID, __VA_ARGS__ };                                              \
                                                                                      \
    /** @brief Map from string names to policy enum values */                         \
    static const std::map<std::string, _type> _dict = {KV_TERMS(__VA_ARGS__)};        \
                                                                                      \
    /** @brief Prints available policy values to stdout */                            \
    static inline void _help() {                                                      \
        std::cout << "Available values for " << #Policy << ":\n";                     \
        for (const auto& entry : _dict) { std::cout << "  " << entry.first << "\n"; } \
    }                                                                                 \
                                                                                      \
    /** @brief Converts string to policy enum value                                   \
     *  @param s String to convert                                                    \
     *  @return Corresponding enum value, or INVALID if not found                     \
     */                                                                               \
    static inline _type _from(const std::string& s) {                                 \
        try {                                                                         \
            auto it = _dict.find(s);                                                  \
            if (it != _dict.end()) { return it->second; }                             \
        } catch (const std::out_of_range& e) { /* Fall through to INVALID */          \
        }                                                                             \
        _help(); /* Print available values for user feedback */                       \
        return INVALID;                                                               \
    }                                                                                 \
    };  // namespace Policy

#endif  // PSND_POLICY_H

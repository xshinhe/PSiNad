#ifndef POLICY_H
#define POLICY_H

#include <map>
#include <string>

#define VAR_NAME(x) #x
#define CAT_NAME(x, y) x##y
#define SELECT_NAME(NAME, NUM) CAT_NAME(NAME##_, NUM)
#define ARG_COUNT(...)                                                  \
    ARG_COUNT_PRIVATE_IMPL(0, ##__VA_ARGS__,                       /**/ \
                           50, 49, 48, 47, 46, 45, 44, 43, 42, 41, /**/ \
                           40, 39, 38, 37, 36, 35, 34, 33, 32, 31, /**/ \
                           30, 29, 28, 27, 26, 25, 24, 23, 22, 21, /**/ \
                           20, 19, 18, 17, 16, 15, 14, 13, 12, 11, /**/ \
                           10, 9, 8, 7, 6, 5, 4, 3, 2, 1,          /**/ \
                           0)

#define ARG_COUNT_PRIVATE_IMPL(_0,                                               /**/ \
                               _1, _2, _3, _4, _5, _6, _7, _8, _9, _10,          /**/ \
                               _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, /**/ \
                               _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, /**/ \
                               _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, /**/ \
                               _41, _42, _43, _44, _45, _46, _47, _48, _49, _50, /**/ \
                               count, ...)                                            \
    count
#define VA_SELECT(NAME, ...) SELECT_NAME(NAME, ARG_COUNT(__VA_ARGS__))(__VA_ARGS__)

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
#define KV_TERMS_1(TERM1, ...) \
    { VAR_NAME(TERM1), TERM1 }

//     , KV_TERMS_0(__VA_ARGS__)
// #define KV_TERMS_0(TERM1, ...) \
//     { VAR_NAME(TERM1), TERM1 }


#define KV_TERMS(...) VA_SELECT(KV_TERMS, ##__VA_ARGS__)

#define DEFINE_POLICY(Policy, ...)                                             \
    namespace Policy {                                                         \
    enum _type { __VA_ARGS__ };                                                \
    static const std::map<std::string, _type> _dict = {KV_TERMS(__VA_ARGS__)}; \
    static inline void _help() {                                               \
        std::cout << "Helps for " << #Policy << ":\n";                         \
        for (auto& i : _dict) std::cout << i.first << " [available]\n";        \
    }                                                                          \
    static inline _type _from(std::string s) {                                 \
        try {                                                                  \
            return _dict.at(s);                                                \
        } catch (std::out_of_range & e) { _help(); }                           \
        return _type(0);                                                       \
    }                                                                          \
    };  // namespace Policy

#endif  // POLICY_H

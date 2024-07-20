#ifndef KIDS_INTRINSIC_RULES_H
#define KIDS_INTRINSIC_RULES_H

namespace PROJECT_NS {

constexpr std::string KENV_OCCUPATION = "$SLR";
constexpr std::string KENV_OCCUPATION = "$MPI";
constexpr std::string KENV_OCCUPATION = "$TRJ";
constexpr std::string KENV_OCCUPATION = "$OCC";

// define some placeholders, where it can read from ENV to get the values

constexpr std::string OCCUPATION_PLACEHOLDER = "$OCC";
constexpr std::string TRJ_INDEX_PLACEHOLDER  = "$TRJ";
constexpr std::string MPI_INDEX_PLACEHOLDER  = "$MPI";

// define some rules, for population transfer, linear spectrum and reaction rate

constexpr std::string KRULE_KKCMM = "{\"rule\":\"KKCMM<mik>(K1<m[$OCC][$OCC]>)\"}";
constexpr std::string KRULE_KKwMM = {"$KK_wMM", "KK"};
constexpr std::string KRULE_KKSQC = {"$KK_SQC", "KK"};

};  // namespace PROJECT_NS

#endif  // KIDS_INTRINSIC_RULES_H
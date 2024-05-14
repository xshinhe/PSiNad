#ifndef KIDS_INTRINSIC_RULES_H
#define KIDS_INTRINSIC_RULES_H

namespace PROJECT_NS {

// reserved notation:

constexpr std::string OCCUPATION_PLACEHOLDER = "$OCC";
constexpr std::string TRJ_INDEX_PLACEHOLDER  = "$TRJ";
constexpr std::string MPI_INDEX_PLACEHOLDER  = "$MPI";

constexpr std::string KRULE_KKCMM = "{\"rule\":\"KKCMM<mik>(K1<m[$OCC][$OCC]>)\"}";
constexpr std::string KRULE_KKwMM = {"$KK_wMM", "KK"};
constexpr std::string KRULE_KKSQC = {"$KK_SQC", "KK"};

};  // namespace PROJECT_NS

#endif  // KIDS_INTRINSIC_RULES_H
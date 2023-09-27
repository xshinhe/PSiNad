#ifndef FLAGS_DICT_H
#define FLAGS_DICT_H

#include <map>

namespace elec_eom {
enum _enum { eac, eaccv, eacfb, rho };
const std::map<std::string, _enum> _dict = {
    {"#eac", eac},      // electronic amplititude c(t)
    {"#eacfb", eacfb},  // forward and backward cf(t), cb(t)
    {"#eaccv", eaccv},  // c(t) and commutator variables gmat(t)
    {"#rho", rho},      // rho(t), @NOTE: all above has relation to rho(t)
};
};  // namespace elec_eom


namespace nad_mem {
enum _enum { def, less, lazy, full };
const std::map<std::string, _enum> _dict = {
    {"#def", def},    // correlation of occupied stated with population array
    {"#less", less},  // correlation of occupied stated with rho matrix
    {"#lazy", lazy},  // correlation of rho matrix with rho matrix
    {"#full", full},  // correlation of rho matrix with rho matrix
};
};  // namespace nad_mem

namespace nad_tcf {
enum _enum { pop, rho, all, tcf, redtcf };
const std::map<std::string, _enum> _dict = {
    {"#pop", pop},       // correlation of occupied stated with population array
    {"#rho", rho},       // correlation of occupied stated with rho matrix
    {"#all", all},       // correlation of rho matrix with rho matrix
    {"tcf", tcf},        // correlation of rho matrix with rho matrix
    {"redtcf", redtcf},  // correlation of rho matrix with rho matrix
};
};  // namespace nad_tcf

namespace elec_init {
enum _enum { occ, eac, rho, d2a, eig, dia, rot, def };
const std::map<std::string, _enum> _dict = {
    {"#occ", occ},  //
    {"#eac", eac},  //
    {"#rho", rho},  //
    {"#d2a", d2a},  //
    {"#eig", eig},  //
    {"#rot", rot},  //
    {"#def", def}   //
};
};  // namespace elec_init


#endif  // FLAGS_DICT_H
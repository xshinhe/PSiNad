/**
 * @file param_wrapper.h
 * @author xshinhe
 * @version 1.0
 * @date 2022-04
 * @brief tools for handle parameter input and access
 * @details
 *      unified interface of parameter struct to different file format. This file
 *      realizes interfaces:
 *
 *  1) types:
 *      * ::Param
 *          * it will be passed into Model class and Solver class.
 *          * it should be accessed to its sub-nodes with indexing method `[]`.
 *      * ::Param_Exception
 *          * it helps to handle errors.
 *
 *  2) methods (or macros):
 *      * Param_ParseFile(): parse file to Param type
 *      * Param_IfHaveKey(): ask key in Param
 *      * Param_FindByKey(): find values in Param type (pre-types and types)
 *      * operator <<()    : overloaded for Param
 *
 *  now support file format:
 *      - [o] JSON (based configor    : https://github.com/Nomango/configor)
 *      - [x] YAML (based yamlcpp     : https://github.com/jbeder/yaml-cpp
 *                       or rapidyaml   : https://github.com/biojppm/rapidyaml)
 *      - [x] TOML (based toml11      : https://github.com/ToruNiina/toml11)
 *      - [x] XML  (based srlzio      : https://github.com/Mazrog/srlzio)
 *
 *  consider later (ENO | namelist) format.
 *
 *  @note
 *
 *      * variable type is labeled [type]
 *      * variable name is labeled <name>
 *      * variable value is labeled {data}
 *      * variable key is labeled /key/
 *
 *      (JSON/YAML/TOML/XML) format will build thier type from reading raw data (we called pre-type of the key),
 *      and we need mapping it to the true-type of the variable in the code. This is what this wrapper does.
 *
 */

#ifndef PARAM_WRAPPER_H
#define PARAM_WRAPPER_H

#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "definitions.h"
#include "phys.h"

#define VAR_NAME(x)            #x
#define CAT_NAME(x, y)         x##y
#define SELECT_NAME(NAME, NUM) CAT_NAME(NAME##_, NUM)
#define ARG_COUNT(...) \
    ARG_COUNT_PRIVATE_IMPL(0, ##__VA_ARGS__, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#define ARG_COUNT_PRIVATE_IMPL(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, count, ...) \
    count
#define VA_SELECT(NAME, ...) SELECT_NAME(NAME, ARG_COUNT(__VA_ARGS__))(__VA_ARGS__)


#if (!defined(PARAM_USE_FORMAT))
#define PARAM_USE_JSON
#else
#if PARAM_USE_FORMAT == "JSON"
#define PARAM_USE_JSON
#elif PARAM_USE_FORMAT == "YAML"
#define PARAM_USE_YAML
#elif PARAM_USE_FORMAT == "TOML"
#define PARAM_USE_TOML
#elif PARAM_USE_FORMAT == "XML"
#define PARAM_USE_XML
#else
#error "Unknown PARAM_USE_FORMAT"
#endif
#endif  // PARAM_USE_FORMAT

/*===================================================
=            Param in use of JSON format            =
===================================================*/
#ifdef PARAM_USE_JSON


#include "../thirdpart/configor/json.hpp"

/**
 *
 * Following is an example of realization of interface to JSON format. Here the json
 * library is used with https://github.com/Nomango/configor.
 *
 * It needs interfaces of:
 *
 * 1) types:
 *      * Param
 *          * it will be passed into Model class and Solver class.
 *          * it should be accessed to its sub-nodes with indexing method `[]`.
 *      * Param_Exception
 *          * it helps to handle errors.
 *
 * 2) methods (or macros):
 *      * Param_ParseFile(): parse file to Param type
 *      * Param_IfHaveKey(): ask key in Param
 *      * Param_FindByKey(): find values in Param type (pre-types and types)
 *      * operator <<()    : overloaded for Param
 */
typedef configor::json Param;
typedef configor::configor_exception Param_Exception;

inline int Param_ParseFile(Param& PM, const std::string& file) {
    std::ifstream ifs(file);
    if (!ifs.is_open()) LOG(FATAL) << "Cannot open file: " << file;
    try {
        ifs >> PM;
    } catch (const Param_Exception& e) { LOG(FATAL) << "Cannot parse file " << file; }
    ifs.close();
    return 0;
}

inline std::ostream& operator<<(std::ostream& os, Param& PM) {
    os << PM.dump(4, ' ');
    return os;
}

inline const char* const_basename(const char* filepath) {
    const char* base = strrchr(filepath, '/');
#ifdef _WIN32  // windows system
    if (!base) base = strrchr(filepath, '\\');
#endif
    return base ? (base + 1) : filepath;
}

#define LOC() (std::string(const_basename(__FILE__)) + ":" + std::to_string(__LINE__))

template <typename T>
inline std::string TypeFLAG() {
    return "unknow";
}
template <>
inline std::string TypeFLAG<bool>() {
    return "bool";
}
template <>
inline std::string TypeFLAG<int>() {
    return "integr";
}
template <>
inline std::string TypeFLAG<double>() {
    return "double";
}
template <>
inline std::string TypeFLAG<std::string>() {
    return "string";
}


#define Param_Report(STATUS, LOCINFO, TYPE, KEY, DATA, MSG)   \
    LOG(STATUS) << LOCINFO << " "           /**/              \
                << "Type <" << TYPE << "> " /*vairable type*/ \
                << "KEY  /" << KEY << "/ "  /*variable key*/  \
                << "DATA {" << DATA << "} : " << MSG;


inline bool Param_IfHaveKey(const Param& PM, const std::string& key) { return !(PM.count(key) == 0); }


template <typename T>
T Param_FindByKey0(const Param& PM, const std::string& key,
                   const std::string& LOCINFO,                  // location name
                   const bool& with_default,                    // if with default choice
                   const T& qval_default,                       // default value
                   const phys::dimension7& qdim = phys::none_d  // default dimension
) {
    T q;
    if (Param_IfHaveKey(PM, key)) {
        switch (PM[key].type()) {
            case configor::config_value_type::string: {
                if (std::is_same<T, std::string>::value) {
                    q = PM[key].get<T>();
                } else if (std::is_same<T, double>::value) {
                    phys::uval uv = phys::us::parse(PM[key].as_string());
                    double qval   = (qdim == phys::none_d) ? phys::au::as(uv.dim, uv) : phys::au::as(qdim, uv);

                    // [double qval] convert to [T q], convert it by stringstream (stupid)
                    std::stringstream sstr;
                    sstr << std::setiosflags(std::ios::scientific) << std::setprecision(32) << qval;
                    sstr >> q;
                } else {  // ERROR
                    Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM[key], "converting fatal");
                }
                break;
            }
            case configor::config_value_type::boolean: {
                if (std::is_same<T, bool>::value) {
                    q = PM[key].get<T>();
                } else {  // ERROR
                    Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM[key], "converting fatal");
                }
                break;
            }
            case configor::config_value_type::number_float: {
                if (std::is_same<T, double>::value) {
                    q = PM[key].as_float();
                } else {  // ERROR
                    Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM[key], "converting fatal");
                }
                break;
            }
            case configor::config_value_type::number_integer: {
                if (std::is_same<T, double>::value) {
                    q = PM[key].as_float();
                } else if (std::is_same<T, int>::value) {
                    q = PM[key].get<T>();
                } else {
                    Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM[key], "converting fatal");
                }
                break;
            }
            default: {  // ERROR
                Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM[key], "converting fatal");
            }
        }
    } else if (with_default) {
        q = qval_default;
        Param_Report(WARNING, LOCINFO, TypeFLAG<T>(), key, qval_default, "use default");
    } else {  // ERROR
        Param_Report(FATAL, LOCINFO, TypeFLAG<T>(), key, PM, "converting fatal");
    }
    return q;
}

#define Param_FindByKey(...)                    VA_SELECT(Param_FindByKey, __VA_ARGS__)
#define Param_FindByKey_4(T, PM, key, dim)      ({ Param_FindByKey0<T>(PM, key, LOC(), false, T(), dim); })
#define Param_FindByKey_5(T, PM, key, dim, def) ({ Param_FindByKey0<T>(PM, key, LOC(), true, def, dim); })



/*=====  End of Param in use of JSON format  ======*/
#endif  // PARAM_USE_XXXX

//===============================================================================================================
// @interfaces: please use interfaces here !!!

#define Param_GetT(T, PM, key, ...) ({ Param_FindByKey(T, PM, key, phys::none_d, ##__VA_ARGS__); })
#define Param_GetV(V, PM, ...)      ({ V = Param_FindByKey(decltype(V), PM, #V, phys::none_d, ##__VA_ARGS__); })
#define Param_GetQ(Q, PM, key, ...) ({ Param_FindByKey(double, PM, key, Q, ##__VA_ARGS__); })

#define Param_Reset(A, B) ({ A = (decltype(A)) B; })

#endif  // PARAM_WRAPPER_H

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

#define Param_ParseFile(PM, file)                                                                         \
    ({                                                                                                    \
        try {                                                                                             \
            std::ifstream ifs(file);                                                                      \
            if (!ifs.is_open()) LOG(FATAL) << "Cannot open file: " << file;                               \
            ifs >> PM;                                                                                    \
            ifs.close();                                                                                  \
        } catch (const Param_Exception& e) { LOG(FATAL) << "Cannot parse file " << file << " to Param"; } \
    })

inline bool Param_IfHaveKey(const Param& PM, const std::string& key) { return !(PM.count(key) == 0); }

template <typename T>
T Param_FindByKey0(const Param& PM, const std::string& key, const phys::dimensions& qdim = phys::none_d) {
    T q;
    switch (PM[key].type()) {
        case configor::config_value_type::string: {
            if (std::is_same<T, double>::value) {  // string-to-double conversion
                phys::uval uv = phys::us::parse_pair(PM[key].as_string());
                double qval   = (qdim == phys::none_d) ? phys::au::as(uv.dim, uv) : phys::au::as(qdim, uv);
                std::stringstream sstr;
                sstr << std::setiosflags(std::ios::scientific) << std::setprecision(32) << qval;
                sstr >> q;
            } else {  // try to get value or throw exception
                q = PM[key].get<T>();
            }
            break;
        }
        case configor::config_value_type::number_float: {
            q = PM[key].as_float();
            break;
        }
        case configor::config_value_type::number_integer: {
            if (std::is_same<T, double>::value) {
                q = PM[key].as_float();
            } else {
                q = PM[key].get<T>();
            }
            break;
        }
        default: {
            q = PM[key].get<T>();
            break;
        }
    }
    return q;
}

#define Param_Report(STATUS, Tname, Kname, DATA, MSG)          \
    LOG(STATUS) << "Type [" << Tname << "] " /*vairable type*/ \
                << "Key /" << Kname << "/ "  /*variable key*/  \
                << "Data {" << DATA << "} " << MSG;

#define Param_FindByKey(...) VA_SELECT(Param_FindByKey, __VA_ARGS__)

// parse 4 parameters
#define Param_FindByKey_4(T,   /*typename*/                                                         \
                          PM,  /*Param object*/                                                     \
                          key, /*key string*/                                                       \
                          dim  /*dimensions*/                                                       \
)                                                                                                   \
    ({                                                                                              \
        T q;                                                                                        \
        if (Param_IfHaveKey(PM, key)) {                                                             \
            try {                                                                                   \
                q = Param_FindByKey0<T>(PM, key, dim);                                              \
            } catch (Param_Exception & e) { Param_Report(FATAL, #T, key, PM[key], "Parse Error"); } \
        } else {                                                                                    \
            Param_Report(FATAL, #T, key, PM, "Is Not Found");                                       \
        }                                                                                           \
        q;                                                                                          \
    })

// parse 5 parameters
#define Param_FindByKey_5(T,   /*typename*/                                                         \
                          PM,  /*Param object*/                                                     \
                          key, /*key string*/                                                       \
                          dim, /*dimensions*/                                                       \
                          def  /*default value*/                                                    \
)                                                                                                   \
    ({                                                                                              \
        T q;                                                                                        \
        if (Param_IfHaveKey(PM, key)) {                                                             \
            try {                                                                                   \
                q = Param_FindByKey0<T>(PM, key, dim);                                              \
            } catch (Param_Exception & e) { Param_Report(FATAL, #T, key, PM[key], "Parse Error"); } \
        } else {                                                                                    \
            q = def;                                                                                \
            Param_Report(WARNING, #T, key, def, "As Default");                                      \
        }                                                                                           \
        q;                                                                                          \
    })

inline std::ostream& operator<<(std::ostream& os, Param& PM) {
    os << PM.dump(4, ' ');
    return os;
}

/*=====  End of Param in use of JSON format  ======*/
#endif  // PARAM_USE_XXXX



//===============================================================================================================
// @interfaces: please use interfaces here !!!

#define Param_GetT(T, PM, key, ...) ({ Param_FindByKey(T, PM, key, phys::none_d, ##__VA_ARGS__); })
#define Param_GetV(V, PM, ...)      ({ V = Param_FindByKey(decltype(V), PM, #V, phys::none_d, ##__VA_ARGS__); })
#define Param_GetQ(Q, PM, key, ...) ({ Param_FindByKey(double, PM, key, Q, ##__VA_ARGS__); })

#define Param_Reset(A, B)                                                                   \
    ({                                                                                      \
        try {                                                                               \
            A = (decltype(A)) B;                                                            \
            LOG(WARNING) << "param /" #A "/ reset to {" << B << "}";                        \
        } catch (const std::runtime_error& e) { LOG(FATAL) << "Reset Error for <" #A ">"; } \
    })

#endif  // PARAM_WRAPPER_H

/**
 * @file Param.h
 * @author xshinhe
 * @date 2023-04
 * @brief class for parameters from input
 * @details
 *  current Param class is based on json format, built on a thirdpart library:
 *  https://github.com/Nomango/configor
 *  herein the used version is: https://github.com/Nomango/configor/tree/2.x
 *
 * @further may support file format:
 *      - [x] YAML (based yamlcpp     : https://github.com/jbeder/yaml-cpp
 *                       or rapidyaml : https://github.com/biojppm/rapidyaml)
 *      - [x] TOML (based toml11      : https://github.com/ToruNiina/toml11)
 *      - [x] XML  (based srlzio      : https://github.com/Mazrog/srlzio)
 *      - [x] ENO or namelist format
 */

#ifndef PARAM_H
#define PARAM_H

#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "../thirdpart/configor/json.hpp"
#include "Exception.h"
#include "phys.h"

/**
 *
 * Following is an example of realization of interface to JSON format.
 * the json library is used with https://github.com/Nomango/configor.
 * the used version is: https://github.com/Nomango/configor/tree/2.x
 */
typedef configor::json JSON;
typedef configor::configor_exception JSON_Exception;

/**
 * @brief return the base of filename
 */
inline const char *basename(const char *filepath) {
    const char *base = strrchr(filepath, '/');
#ifdef _WIN32  // windows system
    if (!base) base = strrchr(filepath, '\\');
#endif
    return base ? (base + 1) : filepath;
}

/**
 * @brief the macro for the location info
 */
#define LOC() (std::string(basename(__FILE__)) + ":" + std::to_string(__LINE__))

namespace PROJECT_NS {

// helpers for types in json
template <typename T>
inline std::string TypeFLAG();
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

/**
 * @brief class Param is a wrapper based on the thirdpart json library
 */
class Param final {
   public:
    enum LoadOption { fromString, fromFile };

    /**
     * @brief constructor
     *
     * @param input : string or filename
     * @param option : Param::fromFile or Param::fromString
     */
    Param(const std::string &input, LoadOption option) : created{true} {
        pj = new JSON();
        // clang-format off
        switch (option) {
            case fromFile: {
                std::ifstream ifs(input);
                try {
                    ifs >> *pj;
                } catch (const JSON_Exception& e) {
                    std::cout << "Invalid file: " << input;
                }
                ifs.close();
                break;
            }
            case fromString: {
                std::stringstream sstr(input);
                try {
                    sstr >> *pj;
                } catch (const JSON_Exception& e) {
                    std::cout << "Invalid string: " << input;
                }
                break;
            }
        }
        // clang-format on
    }

    /**
     * @brief deconstructor (only from the root)
     */
    virtual ~Param() {
        if (created) delete pj;
    }

    // Param &operator[](const std::string &key) {
    //     if (!has_key(key)) { throw param_illegal_key_error(key); }
    //     Param ref = Param();
    //     ref.pj    = &((*pj)[key]);
    //     return ref;
    // }

    bool has_key(const std::string &key) { return !(pj->count(key) == 0); }

    /**
     * @brief get parameter
     * @param key : key of the parameter
     * @param loc : location tracer (always leaves it `LOC()`)
     * @param qdim : physics dimension for conversion in double in AU unit
     * @param default_value: default value is avialable
     * @return T type
     */
    template <typename T, bool Require = false>
    T get(const std::string &key, const std::string &loc, const phys::dimension7 &qdim, const T &default_value) {
        // cannot find the key
        if (!has_key(key)) {
            if (Require) {
                throw param_illegal_key_error(                   //
                    utils::concat(loc,                           //
                                  " Type<", TypeFLAG<T>(), ">",  //
                                  " Key/", key, "/",             //
                                  " : Illegal default")          //
                );
                return T();
            } else {
                try {
                    throw param_warning(                             //
                        utils::concat(loc,                           //
                                      " Type<", TypeFLAG<T>(), ">",  //
                                      " Key/", key, "/",             //
                                      " : Use default ",             //
                                      default_value)                 //
                    );
                } catch (param_warning &w) { std::cerr << w.what() << "\n"; }
                return default_value;
            }
        }
        // the case find the key
        switch ((*pj)[key].type()) {
            case configor::config_value_type::string: {
                if (std::is_same<T, std::string>::value) {
                    return (*pj)[key].get<T>();
                } else if (std::is_same<T, double>::value) {
                    // parse unit
                    phys::uval uv = phys::us::parse((*pj)[key].as_string());
                    double qval   = phys::au::as(qdim, uv);

                    // conversion by stringstream (stupid)
                    T q;
                    std::stringstream ss;
                    ss << std::setiosflags(std::ios::scientific)  //
                       << std::setprecision(32) << qval;          //
                    ss >> q;
                    return q;
                }
                break;
            }
            case configor::config_value_type::boolean: {
                if (std::is_same<T, bool>::value) return (*pj)[key].get<T>();
                break;
            }
            case configor::config_value_type::number_float: {
                T q;
                if (std::is_same<T, double>::value) {
                    q = (*pj)[key].as_float();
                    return q;
                }
                break;
            }
            case configor::config_value_type::number_integer: {
                T q;
                if (std::is_same<T, double>::value) {
                    q = (*pj)[key].as_float();
                    return q;
                } else if (std::is_same<T, int>::value) {
                    return (*pj)[key].get<T>();
                }
                break;
            }
        }
        // cannot be adapted to existing conversions
        throw param_mismatched_type_error(               //
            utils::concat(loc,                           //
                          " Type<", TypeFLAG<T>(), ">",  //
                          " Key/", key, "/",             //
                          " Data{",                      //
                          (*pj)[key].dump(4, ' '), "}",  //
                          " : Converting fatal")         //
        );
        return T();
    }

    // template <typename T>
    // Param &get(KParam<T> *KP, const std::string &loc, const phys::dimension7 &qdim, const T &default_value) {
    //     KP->ptr = get<T>(KP->key, loc, qdim, default_value);
    //     return *this;
    // }

    /**
     * @brief simplified interface
     */
    template <typename T, bool Require = true>
    T get(const std::string &key, const std::string &loc, const phys::dimension7 &qdim) {
        return get<T, Require>(key, loc, qdim, T());
    }
    template <typename T, bool Require = false>
    T get(const std::string &key, const std::string &loc, const T &default_value) {
        return get<T, Require>(key, loc, phys::none_d, default_value);
    }
    template <typename T, bool Require = true>
    T get(const std::string &key, const std::string &loc) {
        return get<T, Require>(key, loc, phys::none_d, T());
    }
    template <typename T, bool Require = true>
    T get(const std::string &key) {
        return get<T, Require>(key, "__loc__", phys::none_d, T());
    }

    /**
     * @pity: Maybe it is a better manner but only suits std=c++2a:
     *
     *  template <typename T>
     *  double get(const std::string& key, const std::string& loc = "", const
     *  T& default_value = T());
     *
     *  template <phys::dimension7 qdim>
     *  double get(const std::string& key, const std::string& loc = "", const
     *  double& default_value = double());
     *
     */


    inline JSON *pjson() { return pj; }

    inline std::string repr() { return pj->dump(4, ' '); }

   private:
    bool created = false;
    JSON *pj;

    Param() : created{false} {};
};

};  // namespace PROJECT_NS

#endif  // PARAM_H

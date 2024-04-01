/**@file        Param.h
 * @brief       Provide struct and interfaces for input parameters
 * @details
 *
 * ## Easier interfaces for getting parameters
 * @note it is a better manner but only suits std=c++2a:
 * ```cpp
 *  template <typename T>
 *  double get(const std::string& key, const std::string& loc = "", const
 *  T& default_value = T());
 *
 *  template <phys::dimension7 qdim>
 *  double get(const std::string& key, const std::string& loc = "", const
 *  double& default_value = double());
 * ```
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> Initial version. Added detailed commentary by ChatGPT.
 * </table>
 **********************************************************************************
 */

#ifndef KIDS_PARAM_H
#define KIDS_PARAM_H

#include <cstring>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#include "../thirdpart/configor/json.hpp"
#include "Exception.h"
#include "phys.h"

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
 * show the location information for debug
 */
#define LOC() (std::string(basename(__FILE__)) + ":" + std::to_string(__LINE__))

namespace PROJECT_NS {

/**
 * this class provides an interface wrapper for the parameter data.
 *
 * @note As a realization of interface to JSON format.
 * the json library is used with https://github.com/Nomango/configor.
 * the used version is: https://github.com/Nomango/configor/tree/2.x
 *
 * @todo More formats should be supported in future, including:
 * - YAML (based yamlcpp     : https://github.com/jbeder/yaml-cpp
 *       or rapidyaml : https://github.com/biojppm/rapidyaml)
 * - TOML (based toml11      : https://github.com/ToruNiina/toml11)
 * - XML  (based srlzio      : https://github.com/Mazrog/srlzio)
 * - ENO or namelist (in Fortran)
 */
class Param final {
   public:
    using JSON           = configor::json;
    using JSON_Exception = configor::configor_exception;

    /** @see Param(const std::string &input, LoadOption option)
     */
    enum LoadOption {
        fromString,  ///< construct Param from string
        fromFile     ///< construct Param from file
    };

    /**
     * @param[in]  input              a Json-like string or a file path
     * @param[in]  option             enum for parse input
     * @ref fromFile                  parse Param from file    \n
     * @ref fromString                parse Param from string  \n
     *
     * ### Usages
     *
     * - initialize Param from string
     *   ```cpp
     *   Param("{\"key\": 1234}", Param::fromString);
     *   ```
     * - initialize Param from file
     *   ```cpp
     *   Param("param_example.json", Param::fromFile);
     *   ```
     */
    Param(const std::string &input, LoadOption option) {
        pj = std::shared_ptr<JSON>(new JSON());
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
     *Param &operator[](const std::string &key) {
     *        if (!has_key(key)) { throw param_illegal_key_error(key); }
     *        Param ref = Param();
     *        ref.pj    = &((*pj)[key]);
     *        return ref;
     *    }
     */

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
                throw kids_error(                              //
                    utils::concat(loc,                         //
                                  " Type<", as_str<T>(), ">",  //
                                  " Key/", key, "/",           //
                                  " : Illegal default")        //
                );
                return T();
            } else {
                try {
                    throw param_warning(                           //
                        utils::concat(loc,                         //
                                      " Type<", as_str<T>(), ">",  //
                                      " Key/", key, "/",           //
                                      " : Use default ",           //
                                      default_value)               //
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
                    phys::uval uv   = phys::us::parse((*pj)[key].as_string());
                    double     qval = phys::au::as(qdim, uv);

                    // conversion by stringstream (stupid)
                    T                 q;
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
        throw kids_error(                                //
            utils::concat(loc,                           //
                          " Type<", as_str<T>(), ">",    //
                          " Key/", key, "/",             //
                          " Data{",                      //
                          (*pj)[key].dump(4, ' '), "}",  //
                          " : Converting fatal")         //
        );
        return T();
    }

    /** @deprecated
     * ```
     * template <typename T>
     * Param &get(KParam<T> *KP, const std::string &loc, const phys::dimension7 &qdim, const T &default_value) {
     *     KP->ptr = get<T>(KP->key, loc, qdim, default_value);
     *     return *this;
     * }
     * ```
     */

    /// @{
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
    /// @}

    inline std::shared_ptr<JSON> pjson() { return pj; }

    inline std::string repr() { return pj->dump(4, ' '); }

   private:
    std::shared_ptr<JSON> pj;
};

};  // namespace PROJECT_NS

#endif  // KIDS_PARAM_H

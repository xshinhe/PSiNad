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

#include "configor/json.hpp"
#include "kids/phys.h"

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
    Param(const std::string &input, LoadOption option);

    /**
     *Param &operator[](const std::string &key) {
     *        if (!has_key(key)) { throw param_illegal_key_error(key); }
     *        Param ref = Param();
     *        ref.pj    = &((*pj)[key]);
     *        return ref;
     *    }
     */

    /**
     * @param[in]  key   The key
     *
     * @return     True if key, False otherwise.
     */
    bool has_key(const std::string &key);

    std::shared_ptr<JSON> pjson();

    std::string repr();

    bool get_bool(const std::string &key, const std::string &loc, const bool &default_value);
    bool get_bool(const std::string &key, const std::string &loc = "__loc__");

    int get_int(const std::string &key, const std::string &loc, const int &default_value);
    int get_int(const std::string &key, const std::string &loc = "__loc__");

    std::string get_string(const std::string &key, const std::string &loc, const std::string &default_value);
    std::string get_string(const std::string &key, const std::string &loc = "__loc__");

    double get_double(const std::string &key, const std::string &loc, const phys::dimension7 &qdim,
                      const double &default_value = double());
    double get_double(const std::string &key, const std::string &loc, const double &default_value);
    double get_double(const std::string &key, const std::string &loc = "__loc__");

   private:
    std::shared_ptr<JSON> pj;

    /**
     * @brief get parameter
     * @param key : key of the parameter
     * @param loc : location tracer (always leaves it `LOC()`)
     * @param qdim : physics dimension for conversion in double in AU unit
     * @param default_value: default value is avialable
     * @return T type
     */
    template <typename T, bool Require = false>
    T get(const std::string &key, const std::string &loc, const phys::dimension7 &qdim, const T &default_value);
};

};  // namespace PROJECT_NS

#endif  // KIDS_PARAM_H

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
 *  This software is a product of Xin's PhD research conducted by Professor
 *Liu's Group at the College of Chemistry and Molecular Engineering, Peking
 *University. All rights are reserved by Peking University. You should have
 *received a copy of the GNU Lesser General Public License along with this
 *software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> Initial version. Added detailed commentary by
 *ChatGPT.
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
#include <vector>

#include "kids/Types.h"
#include "kids/fmt.h"
#include "kids/phys.h"

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
    class DataObject {
       public:
        virtual ~DataObject(){};
    };

    /** @see Param(const std::string &input, LoadOption option)
     */
    enum LoadOption {
        fromString,  ///< construct Param from string
        fromFile     ///< construct Param from file
    };
    enum ImplType { JSON, TOML };

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
     * @param[in]  key   The key
     *
     * @return     True if key, False otherwise.
     */
    bool has_key(const std::string &key);

    bool is_object(const std::string &key);

    bool is_array(const std::string &key);

    bool is_bool(const std::string &key);

    bool is_int(const std::string &key);

    bool is_real(const std::string &key);

    bool is_string(const std::string &key);

    void set_bool(const std::string &key, bool val);

    void set_int(const std::string &key, int val);

    void set_real(const std::string &key, kids_real val);

    void set_string(const std::string &key, std::string val);

    void set_bool_ifndef(const std::string &key, bool val);

    void set_int_ifndef(const std::string &key, int val);

    void set_real_ifndef(const std::string &key, kids_real val);

    void set_string_ifndef(const std::string &key, std::string val);


    kids_bool   get_bool(const std::vector<std::string> &keys, const std::string &loc, const kids_bool &default_value);
    kids_bool   get_bool(const std::vector<std::string> &keys, const std::string &loc = "__loc__");
    kids_int    get_int(const std::vector<std::string> &keys, const std::string &loc, const kids_int &default_value);
    kids_int    get_int(const std::vector<std::string> &keys, const std::string &loc = "__loc__");
    std::string get_string(const std::vector<std::string> &keys, const std::string &loc,
                           const std::string &default_value);
    std::string get_string(const std::vector<std::string> &keys, const std::string &loc = "__loc__");
    kids_real   get_real(const std::vector<std::string> &keys, const std::string &loc, const phys::dimension7 &qdim,
                         const kids_real &default_value = kids_real());
    kids_real   get_real(const std::vector<std::string> &keys, const std::string &loc, const kids_real &default_value);
    kids_real   get_real(const std::vector<std::string> &keys, const std::string &loc = "__loc__");

    std::string repr();

   private:
    ImplType                    impl_t;
    std::shared_ptr<DataObject> obj;
};

};  // namespace PROJECT_NS

#endif  // KIDS_PARAM_H

/**@file        types.h
 * @brief       definition of types in the project and some utiles for types
 * @details
 *              - definition of types
 *              - meta-data of types
 *              - casting/conversion utils
 *
 * ```
 * @author      [author]
 * @date        [latest-date]
 * @version     [version]
 * @copyright   [copyright]
 **********************************************************************************
 * @par revision [logs]:
 * <table>
 * <tr><th> Date    <th> Version    <th> Author    <th> Description
 * <tr><td>[date]   <td>[version]   <td>[author]   <td> [commit]
 * </table>
 *
 **********************************************************************************
 */

// #include "generate/version.h"

#ifndef NUM_TYPES_H
#define NUM_TYPES_H

#include <complex> /* using C++'s std::complex other than C */

namespace PROJECT_NS {

class Param;
class DataSet;

enum kids_dtype {
    kids_void_type,     //
    kids_bool_type,     //
    kids_int_type,      //
    kids_real_type,     //
    kids_complex_type,  //
    kids_str_type,      //
    kids_param_type,    //
    kids_dataset_type,  //
};

using kids_void    = void;                  ///< bool type
using kids_bool    = bool;                  ///< bool type
using kids_int     = int;                   ///< precisions for int type
using kids_real    = double;                ///< precisions for real type
using kids_complex = std::complex<double>;  ///< precisions for complex type
using kids_str     = std::string;           ///< string type
using kids_param   = Param;                 // to be defined
using kids_dataset = DataSet;               // to be defined

/**@name types_cast
 * utils for casting types
 */
///@{
template <class Tto, class Tfrom = void>
inline Tto cast(Tfrom value) {
    return static_cast<Tto>(value);
}

template <class Tto>
inline Tto cast(kids_complex value) {
    return static_cast<Tto>(real(value));
}

template <class Tto, class Tfrom>
inline Tto cast_at(void* data, int index = 0) {
    return cast<Tto>(*((Tfrom*) data + index));
}
///@}

/**
 */
template <typename T>
kids_dtype as_enum();

template <typename T>
std::string as_str();

#define DEFINE_BIND_UTILS_FOR(T)     \
    template <>                      \
    inline kids_dtype as_enum<T>() { \
        return T##_type;             \
    }                                \
    template <>                      \
    inline std::string as_str<T>() { \
        return #T;                   \
    }

DEFINE_BIND_UTILS_FOR(kids_void);
DEFINE_BIND_UTILS_FOR(kids_int);
DEFINE_BIND_UTILS_FOR(kids_bool);
DEFINE_BIND_UTILS_FOR(kids_real);
DEFINE_BIND_UTILS_FOR(kids_complex);
DEFINE_BIND_UTILS_FOR(kids_str);
DEFINE_BIND_UTILS_FOR(kids_param);
DEFINE_BIND_UTILS_FOR(kids_dataset);

};  // namespace PROJECT_NS

#endif  // NUM_TYPES_H

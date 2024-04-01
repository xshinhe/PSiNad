/**@file        Types.h
 * @brief       definition of types in the project and some utiles for types
 * @details
 *              - definition of types
 *              - meta-data of types
 *              - casting/conversion utils
 *
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

#ifndef KIDS_TYPES_H
#define KIDS_TYPES_H

#include <complex>
// #include "generate/version.h"

namespace PROJECT_NS {

// Forward declarations
class Param;
class DataSet;

// Enumeration for data types
enum kids_dtype {
    kids_void_type,     ///< Represents void type
    kids_bool_type,     ///< Represents bool type
    kids_int_type,      ///< Represents integer type
    kids_real_type,     ///< Represents real number type
    kids_complex_type,  ///< Represents complex number type
    kids_str_type,      ///< Represents string type
    kids_param_type,    ///< Represents Param type
    kids_dataset_type   ///< Represents DataSet type
};

// Alias declarations for data types
using kids_void    = void;                  ///< Alias for void type
using kids_bool    = bool;                  ///< Alias for bool type
using kids_int     = int;                   ///< Alias for integer type
using kids_real    = double;                ///< Alias for real number type
using kids_complex = std::complex<double>;  ///< Alias for complex number type
using kids_str     = std::string;           ///< Alias for string type
using kids_param   = Param;                 ///< Alias for Param type (to be defined later)
using kids_dataset = DataSet;               ///< Alias for DataSet type (to be defined later)

/**@name types_cast
 * Utility functions for type casting
 */
///@{
// Casts a value from one type to another
template <class Tto, class Tfrom = void>
inline Tto cast(Tfrom value) {
    return static_cast<Tto>(value);
}

// Specialization for casting complex numbers to another type
template <class Tto>
inline Tto cast(kids_complex value) {
    return static_cast<Tto>(real(value));
}

// Casts a value at a specific memory location to another type
template <class Tto, class Tfrom>
inline Tto cast_at(void* data, int index = 0) {
    return cast<Tto>(*((Tfrom*) data + index));
}
///@}

/**
 * Converts a C++ type to its corresponding kids_dtype enumeration.
 */
template <typename T>
kids_dtype as_enum();

/**
 * Converts a C++ type to its string representation.
 */
template <typename T>
std::string as_str();

// Macro for defining utility functions for type binding
#define DEFINE_BIND_UTILS_FOR(T)     \
    template <>                      \
    inline kids_dtype as_enum<T>() { \
        return T##_type;             \
    }                                \
    template <>                      \
    inline std::string as_str<T>() { \
        return #T;                   \
    }

// Define utility functions for each data type
DEFINE_BIND_UTILS_FOR(kids_void);
DEFINE_BIND_UTILS_FOR(kids_int);
DEFINE_BIND_UTILS_FOR(kids_bool);
DEFINE_BIND_UTILS_FOR(kids_real);
DEFINE_BIND_UTILS_FOR(kids_complex);
DEFINE_BIND_UTILS_FOR(kids_str);
DEFINE_BIND_UTILS_FOR(kids_param);
DEFINE_BIND_UTILS_FOR(kids_dataset);

};  // namespace PROJECT_NS

#endif  // KIDS_TYPES_H

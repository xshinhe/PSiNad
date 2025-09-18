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

#ifndef PSND_TYPES_H
#define PSND_TYPES_H

#include <complex>
// #include "generate/version.h"

namespace PROJECT_NS {

// Forward declarations
class Param;
class Node;
class DataSet;

// Enumeration for data types
enum psnd_dtype {
    psnd_void_type,     ///< Represents void type
    psnd_bool_type,     ///< Represents bool type
    psnd_int_type,      ///< Represents integer type
    psnd_real_type,     ///< Represents real number type
    psnd_complex_type,  ///< Represents complex number type
    psnd_str_type,      ///< Represents string type
    psnd_param_type,    ///< Represents Param type
    psnd_dataset_type   ///< Represents DataSet type
};

// Type aliases for data types
using psnd_void    = void;                  ///< Alias for void type
using psnd_bool    = bool;                  ///< Alias for bool type
using psnd_bint    = int;                   ///< Alias for bool2 type
using psnd_int     = int;                   ///< Alias for integer type
using psnd_real    = double;                ///< Alias for real number type
using psnd_complex = std::complex<double>;  ///< Alias for complex number type
using psnd_str     = std::string;           ///< Alias for string type
using psnd_param   = Param;                 ///< Alias for Param type (to be defined later)
using psnd_dataset = DataSet;               ///< Alias for DataSet type (to be defined later)

/**@name types_cast
 * Utility functions for type casting
 */
///@{
/**
 * @brief Casts a value from one type to another.
 *
 * @tparam Tto The target type.
 * @tparam Tfrom The source type. Defaults to void.
 * @param value The value to cast.
 * @return The casted value.
 */
template <class Tto, class Tfrom = void>
inline Tto cast(Tfrom value) {
    return static_cast<Tto>(value);
}

/**
 * @brief Casts a complex number to a real number type.
 *
 * @tparam Tto The target type.
 * @param value The complex number to cast.
 * @return The real part of the complex number casted to the target type.
 */
template <class Tto>
inline Tto cast_from_complex(psnd_complex value) {
    return static_cast<Tto>(std::real(value));
}

template <>
inline psnd_complex cast_from_complex(psnd_complex value) {
    return value;
}

/**
 * @brief Specialization for casting complex numbers to another type.
 *
 * @tparam Tto The target type.
 * @param value The complex number to cast.
 * @return If the target type is complex, returns the original value; otherwise,
 * returns the real part.
 */
template <class Tto>
inline Tto cast(psnd_complex value) {
    return cast_from_complex<Tto>(value);
}

/**
 * @brief Casts a value at a specific memory location to another type.
 *
 * @tparam Tto The target type.
 * @tparam Tfrom The source type.
 * @param data Pointer to the data.
 * @param index The index (default is 0).
 * @return The casted value.
 */
template <class Tto, class Tfrom>
inline Tto cast_at(void* data, int index = 0) {
    return cast<Tto>(*((Tfrom*) data + index));
}
///@}

/**
 * @brief Converts a C++ type to its corresponding psnd_dtype enumeration.
 *
 * @tparam T The type to convert.
 * @return The corresponding psnd_dtype enumeration value.
 */
template <typename T>
inline psnd_dtype as_enum() {
    static_assert(std::is_same<T, psnd_void>::value || std::is_same<T, psnd_bool>::value ||
                      std::is_same<T, psnd_int>::value || std::is_same<T, psnd_real>::value ||
                      std::is_same<T, psnd_complex>::value || std::is_same<T, psnd_str>::value ||
                      std::is_same<T, psnd_param>::value || std::is_same<T, psnd_dataset>::value,
                  "Unsupported type for as_enum");
    if (std::is_same<T, psnd_void>::value)
        return psnd_dtype::psnd_void_type;
    else if (std::is_same<T, psnd_bool>::value)
        return psnd_dtype::psnd_bool_type;
    else if (std::is_same<T, psnd_int>::value)
        return psnd_dtype::psnd_int_type;
    else if (std::is_same<T, psnd_real>::value)
        return psnd_dtype::psnd_real_type;
    else if (std::is_same<T, psnd_complex>::value)
        return psnd_dtype::psnd_complex_type;
    else if (std::is_same<T, psnd_str>::value)
        return psnd_dtype::psnd_str_type;
    else if (std::is_same<T, psnd_param>::value)
        return psnd_dtype::psnd_param_type;
    else /* if (std::is_same<T, psnd_dataset>::value) */
        return psnd_dtype::psnd_dataset_type;
}

/**
 * @brief Converts a C++ type to its string representation.
 *
 * @tparam T The type to convert.
 * @return The string representation of the type.
 */
template <typename T>
inline std::string as_str() {
    if (std::is_same<T, psnd_void>::value)
        return "psnd_void";
    else if (std::is_same<T, psnd_bool>::value)
        return "psnd_bool";
    else if (std::is_same<T, psnd_int>::value)
        return "psnd_int";
    else if (std::is_same<T, psnd_real>::value)
        return "psnd_real";
    else if (std::is_same<T, psnd_complex>::value)
        return "psnd_complex";
    else if (std::is_same<T, psnd_str>::value)
        return "psnd_str";
    else if (std::is_same<T, psnd_param>::value)
        return "psnd_param";
    else if (std::is_same<T, psnd_dataset>::value)
        return "psnd_dataset";
    else
        return "unknown";
}

/**
 * @brief Converts a C++ enum to its string representation.
 *
 * @tparam T The type to convert.
 * @return The string representation of the type.
 */
inline std::string enum_t_as_str(psnd_dtype type) {
    if (type == psnd_dtype::psnd_void_type)
        return "psnd_void";
    else if (type == psnd_dtype::psnd_bool_type)
        return "psnd_bool";
    else if (type == psnd_dtype::psnd_int_type)
        return "psnd_int";
    else if (type == psnd_dtype::psnd_real_type)
        return "psnd_real";
    else if (type == psnd_dtype::psnd_complex_type)
        return "psnd_complex";
    else if (type == psnd_dtype::psnd_str_type)
        return "psnd_str";
    else if (type == psnd_dtype::psnd_param_type)
        return "psnd_param";
    else if (type == psnd_dtype::psnd_dataset_type)
        return "psnd_dataset";
    else
        return "unknown";
}

};  // namespace PROJECT_NS

#endif  // PSND_TYPES_H

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

#ifndef KIDS_TYPES_H
#define KIDS_TYPES_H

#include <complex>
// #include "generate/version.h"

namespace PROJECT_NS {

// Forward declarations
class Param;
class Node;
class DataSet;

// Enumeration for data types
enum kids_dtype
{
    kids_void_type,    ///< Represents void type
    kids_bool_type,    ///< Represents bool type
    kids_int_type,     ///< Represents integer type
    kids_real_type,    ///< Represents real number type
    kids_complex_type, ///< Represents complex number type
    kids_str_type,     ///< Represents string type
    kids_param_type,   ///< Represents Param type
    kids_dataset_type  ///< Represents DataSet type
};

// Type aliases for data types
using kids_void = void;                    ///< Alias for void type
using kids_bool = bool;                    ///< Alias for bool type
using kids_bint = int;                     ///< Alias for bool2 type
using kids_int = int;                      ///< Alias for integer type
using kids_real = double;                  ///< Alias for real number type
using kids_complex = std::complex<double>; ///< Alias for complex number type
using kids_str = std::string;              ///< Alias for string type
using kids_param = Param;     ///< Alias for Param type (to be defined later)
using kids_dataset = DataSet; ///< Alias for DataSet type (to be defined later)

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
inline Tto cast(Tfrom value)
{
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
inline Tto cast_from_complex(kids_complex value)
{
    return static_cast<Tto>(std::real(value));
}

template <>
inline kids_complex cast_from_complex(kids_complex value)
{
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
inline Tto cast(kids_complex value)
{
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
inline Tto cast_at(void* data, int index = 0)
{
    return cast<Tto>(*((Tfrom*)data + index));
}
///@}

/**
 * @brief Converts a C++ type to its corresponding kids_dtype enumeration.
 *
 * @tparam T The type to convert.
 * @return The corresponding kids_dtype enumeration value.
 */
template <typename T>
inline kids_dtype as_enum()
{
    static_assert(std::is_same<T, kids_void>::value ||
                      std::is_same<T, kids_bool>::value ||
                      std::is_same<T, kids_int>::value ||
                      std::is_same<T, kids_real>::value ||
                      std::is_same<T, kids_complex>::value ||
                      std::is_same<T, kids_str>::value ||
                      std::is_same<T, kids_param>::value ||
                      std::is_same<T, kids_dataset>::value,
                  "Unsupported type for as_enum");
    if (std::is_same<T, kids_void>::value)
        return kids_dtype::kids_void_type;
    else if (std::is_same<T, kids_bool>::value)
        return kids_dtype::kids_bool_type;
    else if (std::is_same<T, kids_int>::value)
        return kids_dtype::kids_int_type;
    else if (std::is_same<T, kids_real>::value)
        return kids_dtype::kids_real_type;
    else if (std::is_same<T, kids_complex>::value)
        return kids_dtype::kids_complex_type;
    else if (std::is_same<T, kids_str>::value)
        return kids_dtype::kids_str_type;
    else if (std::is_same<T, kids_param>::value)
        return kids_dtype::kids_param_type;
    else /* if (std::is_same<T, kids_dataset>::value) */
        return kids_dtype::kids_dataset_type;
}

/**
 * @brief Converts a C++ type to its string representation.
 *
 * @tparam T The type to convert.
 * @return The string representation of the type.
 */
template <typename T>
inline std::string as_str()
{
    if (std::is_same<T, kids_void>::value)
        return "kids_void";
    else if (std::is_same<T, kids_bool>::value)
        return "kids_bool";
    else if (std::is_same<T, kids_int>::value)
        return "kids_int";
    else if (std::is_same<T, kids_real>::value)
        return "kids_real";
    else if (std::is_same<T, kids_complex>::value)
        return "kids_complex";
    else if (std::is_same<T, kids_str>::value)
        return "kids_str";
    else if (std::is_same<T, kids_param>::value)
        return "kids_param";
    else if (std::is_same<T, kids_dataset>::value)
        return "kids_dataset";
    else
        return "unknown";
}

/**
 * @brief Converts a C++ enum to its string representation.
 *
 * @tparam T The type to convert.
 * @return The string representation of the type.
 */
inline std::string enum_t_as_str(kids_dtype type)
{
    if (type == kids_dtype::kids_void_type)
        return "kids_void";
    else if (type == kids_dtype::kids_bool_type)
        return "kids_bool";
    else if (type == kids_dtype::kids_int_type)
        return "kids_int";
    else if (type == kids_dtype::kids_real_type)
        return "kids_real";
    else if (type == kids_dtype::kids_complex_type)
        return "kids_complex";
    else if (type == kids_dtype::kids_str_type)
        return "kids_str";
    else if (type == kids_dtype::kids_param_type)
        return "kids_param";
    else if (type == kids_dtype::kids_dataset_type)
        return "kids_dataset";
    else
        return "unknown";
}

}; // namespace PROJECT_NS

#endif // KIDS_TYPES_H

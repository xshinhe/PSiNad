/**@file        Variable.h
 * @brief       Defines the Variable class for managing program variables
 * @details
 *  This file declares the `VARIABLE_BASE` abstract base class and `VARIABLE` template class,
 *  which encapsulate variables used in the program. Each variable contains metadata such as
 *  name, data type, shape (dimensions), and documentation string. The base class provides
 *  a unified interface for variable management, while the template class handles type-specific
 *  implementations.
 *
 * @author      Xin He
 * @date        2024-04
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
 * <tr><td> 2024-04-01  <td> Initial version.
 * <tr><td> [YYYY-MM-DD]<td> Added exception handling and improved documentation.
 * </table>
 *
 **********************************************************************************
 */

#ifndef PSND_VARIABLE_H
#define PSND_VARIABLE_H

#include <string>
#include <vector>

#include "psnd/Exception.h"  // Include exception handling utilities
#include "psnd/Shape.h"

namespace PROJECT_NS {

/**
 * @brief Abstract base class for all variables, providing a unified interface
 * @details Serves as a base class for the templated `VARIABLE` class. Maintains a static list
 *          of all variable instances for global management and introspection.
 */
class VARIABLE_BASE {
   public:
    /**
     * @brief Pure virtual method to get the variable name
     * @return Name of the variable (std::string)
     */
    virtual std::string name() const = 0;

    /**
     * @brief Pure virtual method to get the variable documentation
     * @return Documentation string (std::string)
     */
    virtual std::string doc() const = 0;

    /**
     * @brief Static list storing pointers to all variable instances
     * @details Facilitates global enumeration and management of variables.
     */
    static std::vector<VARIABLE_BASE*> _LIST;
};

// Initialize the static list (see in vars_list.cpp)
// inline std::vector<VARIABLE_BASE*> VARIABLE_BASE::_LIST;

/**
 * @brief Helper function to replace substrings in a string
 * @param resource_str Original string to modify
 * @param sub_str Substring to be replaced
 * @param new_str New substring to insert
 * @return Modified string with all occurrences of `sub_str` replaced by `new_str`
 * @warning Throws `psnd_error` if `sub_str` is empty (to avoid infinite loops)
 */
inline std::string _subreplace(std::string resource_str, std::string sub_str, std::string new_str) {
    // Validate input to prevent infinite loops
    if (sub_str.empty()) { throw psnd_error("Invalid substring: cannot replace empty string in _subreplace"); }

    std::string            dst_str = resource_str;
    std::string::size_type pos     = 0;
    while ((pos = dst_str.find(sub_str, pos)) != std::string::npos) {
        dst_str.replace(pos, sub_str.length(), new_str);
        pos += new_str.length();  // Move past the replaced segment to avoid re-matching
    }
    return dst_str;
}

/**
 * @brief Templated class representing a concrete variable with type-specific data
 * @tparam T Data type of the variable (e.g., int, float, std::string)
 * @details Inherits from `VARIABLE_BASE` and implements its interface. Stores variable metadata
 *          (name, shape, documentation) and预留 (reserves) space for type-specific data via `_data`.
 */
template <class T>
class VARIABLE final : public VARIABLE_BASE {
   public:
    /**
     * @brief Constructor to initialize a variable with metadata
     * @param name Name of the variable (will have "::" replaced with "." via `_subreplace`)
     * @param shape Pointer to a `Shape` object defining the variable's dimensions (non-null)
     * @param doc Documentation string describing the variable's purpose
     * @warning Throws `psnd_error` if `shape` is a null pointer (critical for valid dimension checks)
     * @note Automatically adds the instance to `VARIABLE_BASE::_LIST` upon construction
     */
    VARIABLE(const std::string& name, Shape* shape, const std::string& doc)
        : _name{_subreplace(name, "::", ".")}, _shape{shape}, _doc{doc}, allocated(false) {
        // Enforce non-null shape (critical for subsequent shape() calls)
        if (shape == nullptr) { throw psnd_error("VARIABLE constructor: shape cannot be a null pointer"); }
        VARIABLE_BASE::_LIST.push_back(this);
    }

    /**
     * @brief Get the variable's name (with "::" replaced by ".")
     * @return Modified name as a std::string
     */
    std::string name() const override { return _name; }

    /**
     * @brief Get the variable's documentation string
     * @return Documentation as a std::string
     */
    std::string doc() const override { return _doc; }

    /**
     * @brief Get the variable's shape (dimensions)
     * @return Reference to the `Shape` object (guaranteed non-null due to constructor checks)
     */
    Shape& shape() const {
        psnd_assert(_shape != nullptr, "VARIABLE::shape(): _shape is null (invalid state)");
        return (*_shape);
    }

   private:
    std::string _name;      ///< Variable name (with "::" replaced by "." for consistency)
    T*          _data;      ///< Pointer to the variable's data (预留 for future use)
    Shape*      _shape;     ///< Pointer to the `Shape` object defining dimensions (non-null)
    std::string _doc;       ///< Documentation string describing the variable
    bool        allocated;  ///< Flag indicating if `_data` has been dynamically allocated
};

};  // namespace PROJECT_NS

#endif  // PSND_VARIABLE_H

/**@file        Shape.h
 * @brief       Declaration of the Shape class
 * @details     This file provides the declaration of the Shape class. It serves
 *              as contoller for the shape of tensors.
 *
 *              It provides both static shape container std::vector<std::size_t>,
 *              as well as dynamic shape container std::vector<std::size_t*>.
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @todo
 *  support dynamic shape in a more friendly approach (reduce the time cost)
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> initial version. Extract Shape class from the file
 *                          DataSet.h. Make Shape suitable for dynamic dimensions.
 * </table>
 *
 **********************************************************************************
 */

#ifndef PSND_SHAPE_H
#define PSND_SHAPE_H

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

namespace PROJECT_NS {

using Basic_Dimen = std::size_t;               ///< Fundamental type for a single dimension (size_t)
using Basic_Shape = std::vector<std::size_t>;  ///< Type for static shape representation (vector of dimensions)

/**
 * @class Shape
 * @brief Manages tensor shape information with support for static and dynamic dimensions
 * @details
 * The Shape class encapsulates the dimensionality (rank) and size of tensors. It handles two types of shapes:
 * - Static: Dimensions are fixed and stored in a vector of size_t.
 * - Dynamic: Dimensions are tracked via pointers to size_t variables, allowing automatic updates when the pointed
 * values change. It also computes leading dimensions (_ldims) where _ldims[i] = product of dimensions from i+1 to
 * rank-1, which is critical for efficient tensor index calculations.
 */
class Shape {
   public:
    using dynamic_t =
        std::vector<std::size_t*>;  ///< Type for dynamic shape representation (vector of dimension pointers)

    /**
     * @brief Construct a rank-1 static Shape from a single dimension
     * @param size Size of the single dimension (rank = 1)
     */
    Shape(std::size_t size) : _rank{1}, _dims{{size}}, _ldims{1}, _size{size} { enable_dynamic = false; }

    /**
     * @brief Construct a static Shape from a vector of dimensions
     * @param dims Vector of static dimensions (rank = dims.size())
     */
    Shape(std::vector<std::size_t> dims) : _rank{dims.size()}, _dims{dims}, _ldims(dims.size(), 0), _size{1} {
        enable_dynamic = false;
        update();  // Initialize _ldims and _size
    }

    /**
     * @brief Construct a dynamic Shape from pointers to dimensions
     * @param dims_ptr Vector of pointers to size_t variables (dimensions update when pointed values change)
     */
    Shape(std::vector<std::size_t*> dims_ptr)
        : _rank{dims_ptr.size()},
          _dims_ptr{dims_ptr},
          _dims(_dims_ptr.size(), 0),
          _ldims(_dims_ptr.size(), 0),
          _size{1} {
        enable_dynamic = true;
    }

/**
 * @brief Macro to simplify construction of dynamic shapes
 * @details Expands to a Shape instance initialized with a dynamic_t (vector of size_t pointers)
 * @param ... Arguments to construct the dynamic_t vector
 */
#define ShapeDynamic(...) Shape(Shape::dynamic_t{__VA_ARGS__})

    /**
     * @brief Get the rank (number of dimensions) of the shape
     * @details Updates dynamic dimensions (if enabled) before returning the rank
     * @return Rank as an integer (0 for scalar, 1 for vector, etc.)
     */
    inline int rank() {
        if (enable_dynamic) update();  // Sync dynamic dimensions first
        return _rank;
    }

    /**
     * @brief Get the total number of elements described by the shape
     * @details Updates dynamic dimensions (if enabled) before returning the size (product of all dimensions)
     * @return Total number of elements (1 for rank-0, product of dimensions for rank >= 1)
     */
    inline int size() {
        if (enable_dynamic) update();  // Sync dynamic dimensions first
        return _size;
    }

    /**
     * @brief Overload operator() to return the total size
     * @return Same as size() - total number of elements
     */
    inline int operator()() { return size(); }

    /**
     * @brief Get the current dimensions of the shape
     * @details Updates dynamic dimensions (if enabled) before returning the dimension vector
     * @return Reference to vector of current dimensions
     */
    inline std::vector<std::size_t>& dims() {
        if (enable_dynamic) update();  // Sync dynamic dimensions first
        return _dims;
    }

    /**
     * @brief Convert the shape to static mode by freezing current dimensions
     * @details Disables dynamic updates and caches the current dimension values.
     * Subsequent changes to underlying variables (for dynamic shapes) will not affect this Shape.
     */
    inline void static_build() {
        update();                // Final sync of dynamic dimensions
        enable_dynamic = false;  // Disable further updates
    }

    /**
     * @brief Convert the shape to a string representation
     * @return String in the format "<dim1,dim2,...,dimN,>" (e.g., "<2,3,>" for rank-2 shape [2,3])
     */
    inline std::string to_string() {
        if (enable_dynamic) update();  // Sync dynamic dimensions first
        std::stringstream ss;
        ss << "<";
        for (auto& i : _dims) ss << i << ",";
        ss << ">";
        return ss.str();
    }

    /**
     * @brief Check equality with another Shape
     * @return true/false
     */
    inline bool operator==(const Shape& other) const {
        if (enable_dynamic != other.enable_dynamic) return false;
        if (_rank != other._rank) return false;

        // Ensure dynamic shapes are updated before comparison
        if (enable_dynamic) {
            // Create mutable copies to call update()
            Shape temp_this  = *this;
            Shape temp_other = other;
            temp_this.update();
            temp_other.update();
            return temp_this._dims == temp_other._dims && temp_this._size == temp_other._size;
        } else {
            return _dims == other._dims && _size == other._size;
        }
    }

    /**
     * @brief Inequality operator (uses ==)
     * @return true/false
     */
    inline bool operator!=(const Shape& other) const { return !(*this == other); }

   private:
    /**
     * @brief Update dimensions, leading dimensions, and total size
     * @details For dynamic shapes: Syncs _dims with values pointed to by _dims_ptr (validates dimensions > 0).
     * For all shapes: Recalculates leading dimensions (_ldims) and total size (_size).
     * Special handling for rank-0 (scalar) shapes: Sets size to 1.
     */
    void update() {
        if (enable_dynamic) {
            // Sync dynamic dimensions from pointers and validate
            for (int i = 0; i < _rank; ++i) {
                _dims[i] = *(_dims_ptr[i]);  // Dereference pointer to get current dimension
                if (_dims[i] <= 0) throw std::invalid_argument("Dimension cannot be zero or negative");
            }
        }

        // Handle rank-0 (scalar) case to avoid out-of-bounds access
        if (_rank == 0) {
            _size = 1;  // Convention: scalar (rank-0) has size 1
            return;
        }

        // Calculate leading dimensions and total size for rank >= 1
        _ldims[_rank - 1] = 1;                 // Leading dimension of the last rank is always 1
        _size             = _dims[_rank - 1];  // Initialize size with the last dimension

        // Compute leading dimensions for higher ranks (from rank-2 down to 0)
        for (int i = _rank - 2; i >= 0; --i) {
            _ldims[i] = _ldims[i + 1] * _dims[i + 1];  // Product of subsequent dimensions
            _size *= _dims[i];                         // Accumulate total size
        }
    }

    bool                      enable_dynamic = false;  ///< Flag indicating if dynamic dimension updates are enabled
    std::size_t               _rank;                   ///< Number of dimensions (rank) of the shape
    std::vector<std::size_t*> _dims_ptr;  ///< Pointers to dynamic dimensions (used only if enable_dynamic is true)
    std::vector<std::size_t>  _dims;      ///< Current dimensions (synced with _dims_ptr for dynamic shapes)
    std::vector<std::size_t>  _ldims;     ///< Leading dimensions for efficient index calculation
    std::size_t               _size;      ///< Total number of elements (product of all dimensions)
};

};  // namespace PROJECT_NS

#endif  // PSND_SHAPE_H
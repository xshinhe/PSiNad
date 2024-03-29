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
 *  This software is part of the research conducted by the Prof. Liu's Group at the
 *  College of Chemistry and Molecular Engineering (CCME), Peking University.
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

#ifndef KIDS_SHAPE_H
#define KIDS_SHAPE_H

#include <cassert>
#include <vector>

namespace PROJECT_NS {

using Basic_Dimen = std::size_t;
using Basic_Shape = std::vector<std::size_t>;

/**
 * Shape class provide information about of a Tensor's shape
 */
class Shape {
   public:
    ///< construct from a number (rank-1 Shape)
    Shape(std::size_t size) : _rank{1}, _dims{{size}}, _ldims{1}, _size{size} { enable_dynamic = false; }

    ///< construct from a vector (static shape)
    Shape(std::vector<std::size_t> dims) : _rank{dims.size()}, _dims{dims}, _ldims(dims.size(), 0), _size{1} {
        enable_dynamic = false;
        update();
    }

    ///< construct from a vector of std::size_t pointer (dynamic shape)
    Shape(std::vector<std::size_t*> dims_ptr)
        : _rank{dims_ptr.size()},
          _dims_ptr{dims_ptr},
          _dims(_dims_ptr.size(), 0),
          _ldims(_dims_ptr.size(), 0),
          _size{1} {
        enable_dynamic = true;
    }

    ///< get rank of a Shape
    inline int rank() {
        if (enable_dynamic) update();
        return _rank;
    }

    ///< get data size described by a Shape
    inline int size() {
        if (enable_dynamic) update();
        return _size;
    }

    inline std::vector<std::size_t>& dims() {
        if (enable_dynamic) update();
        return _dims;
    }

    inline void static_build() {
        update();
        enable_dynamic = false;
    }

   private:
    void update() {
        if (enable_dynamic)
            for (int i = 0; i < _rank; ++i) {
                _dims[i] = *(_dims_ptr[i]);
                assert(_dims[i] > 0);
            }
        _ldims[_rank - 1] = 1;
        _size             = _dims[_rank - 1];
        for (int i = _rank - 2; i >= 0; --i) {
            _ldims[i] = _ldims[i + 1] * (_dims[i + 1]);
            _size *= _dims[i];
        }
    }

    bool enable_dynamic = false;

    std::size_t               _rank;      ///< rank of the shape
    std::vector<std::size_t*> _dims_ptr;  ///< pointers to dimensions for each rank
    std::vector<std::size_t>  _dims;      ///< dimensions for each rank
    std::vector<std::size_t>  _ldims;     ///< leading dimensions
    std::size_t               _size;      ///< size of data
};

};  // namespace PROJECT_NS

#endif  // KIDS_SHAPE_H

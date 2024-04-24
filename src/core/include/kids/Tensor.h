/**@file        Tensor.h
 * @brief       A class inherited from Node, representing tensors.
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
 *
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-22  <td> move from DataSet.h
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_TENSOR_H
#define KIDS_TENSOR_H

#include <memory>

#include "kids/Node.h"
#include "kids/Shape.h"
#include "kids/Types.h"
#include "kids/fmt.h"

namespace PROJECT_NS {

/**
 * Tensor class is the container for array-like data
 */
template <typename T>
class Tensor final : public Node {
   public:
    using SizeType = std::size_t;
    using DataType = T;

    Tensor(Shape S, const std::string& info = "") : _shape{S}, _doc_info{info} {
        _type = as_enum<T>();
        _size = _shape.size();
        _data = std::shared_ptr<std::vector<T>>(new std::vector<T>(_size, 0));
    }

    virtual std::string repr() {
        std::ostringstream os;
        T*                 ptr = _data->data();
        os << as_str<T>();
        os << FMT(0) << _size;
        os << "\n";
        for (int i = 0; i < _size; ++i) os << FMT(8) << ptr[i];
        return os.str();
    }

    virtual std::string help(const std::string& name) { return _doc_info; }

    inline T* data() { return _data->data(); }

    inline std::size_t size() { return _size; }

    inline Shape& shape() { return _shape; }

   private:
    std::size_t                     _size;
    Shape                           _shape;
    std::string                     _doc_info;
    std::shared_ptr<std::vector<T>> _data;
};

};  // namespace PROJECT_NS

#endif  // KIDS_TENSOR_H
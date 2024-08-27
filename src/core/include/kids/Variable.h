/**@file        Variable.h
 * @brief       this file provide Variable class
 * @details
 *  Variable class defines variables used in the program. It contains type, shape,
 *  and doc_string information.
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
 * <tr><td> 2024-04-01  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Variable_H
#define KIDS_Variable_H

#include <vector>

#include "kids/Shape.h"

namespace PROJECT_NS {

class VARIABLE_BASE {
   public:
    virtual std::string name() const = 0;
    virtual std::string doc() const  = 0;

    static std::vector<VARIABLE_BASE*> _LIST;
};

inline std::string _subreplace(std::string resource_str, std::string sub_str, std::string new_str) {
    std::string            dst_str = resource_str;
    std::string::size_type pos     = 0;
    while ((pos = dst_str.find(sub_str)) != std::string::npos) { dst_str.replace(pos, sub_str.length(), new_str); }
    return dst_str;
}

template <class T>
class VARIABLE final : public VARIABLE_BASE {
   public:
    VARIABLE(const std::string& name, Shape* shape, const std::string& doc)
        : _name{_subreplace(name, "::", ".")}, _shape{shape}, _doc{doc} {
        VARIABLE_BASE::_LIST.push_back(this);
    }

    std::string name() const { return _name; }

    std::string doc() const { return _doc; }

    Shape& shape() const { return (*_shape); }

   private:
    std::string _name;
    T*          _data;
    Shape*      _shape;
    std::string _doc;
    bool        allocated = false;
};

};  // namespace PROJECT_NS

#endif  // KIDS_Variable_H

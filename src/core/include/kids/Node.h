/**@file        Node.h
 * @brief       Declaration of the Node class used for DataSet.
 * @details     This file provides an abstract interface class for tree nodes.
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
 *
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-22  <td> move from DataSet.h
 *
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_NODE_H
#define KIDS_NODE_H

#include "kids/Types.h"

namespace PROJECT_NS {

/**
 * Node class is an abstract interface. A Node will store a tensor or other Nodes
 */
class Node {
   public:
    virtual std::string repr() = 0;

    virtual std::string help(const std::string& name) = 0;

    inline kids_dtype type() { return _type; }

   protected:
    friend class DataSet;

    kids_dtype _type = kids_void_type;
};
};  // namespace PROJECT_NS

#endif  // KIDS_NODE_H

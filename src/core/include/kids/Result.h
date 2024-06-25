/**@file        Result.h
 * @brief       provide Result class
 * @details     for passing selected data in DataSet and RuleSet
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
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-05-14  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_RESULT_H
#define KIDS_RESULT_H

#include <tuple>
#include <vector>

#include "kids/Types.h"

namespace PROJECT_NS {
class Result final {
   public:
    using DataType      = std::vector<std::tuple<std::string,  // data name (without field)
                                            void*,        // data pointer
                                            kids_dtype,   // data type
                                            std::size_t,  // size of data
                                            std::size_t   // stride of data
                                            >>;
    using HeaderType    = std::vector<std::string>;
    using DataframeType = std::vector<std::vector<kids_real>>;

    static Result defineAtNewField(const std::string& field, Result& res_in);

    // std::pair<HeaderType, DataframeType>
    // int asDataframe() { return 0; }
    //

    Result(){};

    DataType& data() { return _data; }

   private:
    friend class RuleSet;
    DataType _data;
};

};  // namespace PROJECT_NS

#endif  // KIDS_RESULT_H

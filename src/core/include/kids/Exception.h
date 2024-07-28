/**@file        Exception.h
 * @brief       provide Exception structs
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
 * <tr><td> 2024-04-01  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef OPENDF_EXCEPTION_H
#define OPENDF_EXCEPTION_H


#include "kids/concat.h"

namespace PROJECT_NS {

const std::string MSG_PARAM_BAD_KEY    = "bad key error : ";
const std::string MSG_PARAM_BAD_TYPE   = "bad type error : ";
const std::string MSG_DATASET_BAD_KEY  = "bad key error : ";
const std::string MSG_DATASET_BAD_SIZE = "bad size error : ";
const std::string MSG_DATASET_BAD_TYPE = "bad type error : ";
const std::string MSG_DATASET_ERR_LOAD = "load error : ";
const std::string MSG_DATASET_ERR_DUMP = "dump error : ";

struct kids_error : public std::runtime_error {
    kids_error(std::string const text) : std::runtime_error(text) {}
};

struct param_warning : public kids_error {
    param_warning(std::string const text) : kids_error(utils::concat("parameter warning of : ", text)) {}
};

#define kids_assert(condition, info) \
    if (!(condition)) throw kids_error(utils::concat(info, "\nassert on: [", #condition, "] fails!"));


};  // namespace PROJECT_NS

#endif  // OPENDF_EXCEPTION_H

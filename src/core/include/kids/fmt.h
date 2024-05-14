/**@file        fmt.h
 * @brief       utils for formating
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
 * <tr><td> 2024-05-01  <td> moved here for initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_FMT_H
#define KIDS_FMT_H

#include <iomanip>
#include <iostream>

namespace PROJECT_NS {

/**
 * control the io printing format
 */
constexpr inline int FMT_WIDTH(int X) { return X + 7; }
#define FMT(X)                                                            \
    " " << std::setiosflags(std::ios::scientific) /*scientific notation*/ \
        << std::setprecision(X)                   /*precision*/           \
        << std::right                             /*alignment*/           \
        << std::setw(FMT_WIDTH(X))                /*width of text*/

/**
 * show the location information for debug
 */
#define LOC() (std::string(basename(__FILE__)) + ":" + std::to_string(__LINE__))

};  // namespace PROJECT_NS

#endif  // KIDS_FMT_H

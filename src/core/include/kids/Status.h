/**@file        Status.h
 * @brief       provide Status class
 * @details     Status
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
 * <tr><td> 2024-05-14  <td> moved here for initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Status_H
#define KIDS_Status_H

namespace PROJECT_NS {

struct Status {
    Status() : succ{true}, stage{0}, mpi_rank{0}, icalc{0} {};

    Status(bool succ, int stage = 0, int mpi_rank = 0, int icalc = 0)
        : succ{succ}, stage{stage}, mpi_rank{mpi_rank}, icalc{icalc} {};

    // Status(const Status&) = delete;  ///< Disable copy constructor

    // Status& operator=(const Status&) = delete;  ///< Disable copy assignment operator

    bool succ     = true;
    int  stage    = 0;
    int  mpi_rank = 0;
    int  icalc    = 0;
};

};  // namespace PROJECT_NS

#endif  // KIDS_Status_H

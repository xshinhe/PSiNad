/**@file        debug_utils.h
 * @brief       provide utils for debugging the code
 * @details     utils printing prettly in the format of numpy array.
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
 * @warning    Do not include this file to any header. You'd better include it only
 *  in source files!
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_DEBUG_UTILS_H
#define KIDS_DEBUG_UTILS_H

#include <iomanip>
#include <iostream>

#include "kids/fmt.h"

#define PRINT_ARRAY(_A, _n1, _n2)                                                           \
    ({                                                                                      \
        std::cout << #_A << " = np.array([\n";                                              \
        int _idxA = 0;                                                                      \
        for (int _i = 0; _i < (_n1); ++_i) {                                                \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++] << ","; \
            std::cout << std::endl;                                                         \
        }                                                                                   \
        { std::cout << "])\n"; }                                                            \
    })

#define PRINT_ARRAY_8(_A, _n1, _n2)                                                         \
    ({                                                                                      \
        std::cout << #_A << " = np.array([\n";                                              \
        int _idxA = 0;                                                                      \
        for (int _i = 0; _i < (_n1); ++_i) {                                                \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(8) << (_A)[_idxA++] << ","; \
            std::cout << std::endl;                                                         \
        }                                                                                   \
        { std::cout << "])\n"; }                                                            \
    })


#endif  // KIDS_DEBUG_UTILS_H

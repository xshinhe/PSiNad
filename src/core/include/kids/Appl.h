/**@file        Appl.h
 * @brief       Declaration of the Appl class
 * @details     This file helps in manage context of Kernels
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
 * <tr><td> 2024-06-23  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Appl_H
#define KIDS_Appl_H

#include "kids/DataSet.h"
#include "kids/RuleSet.h"

namespace PROJECT_NS {

class Appl : public std::enable_shared_from_this<Appl> {
    Appl();
    virtual std::shared_ptr<RuleSet> GenerateRuleSet(std::shared_ptr<DataSet> DS);
};
};  // namespace PROJECT_NS

#endif  // KIDS_Appl_H

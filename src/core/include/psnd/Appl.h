/**@file        Application.h
 * @brief       Declaration of the Application
 * @details     This file defines Applications built over dynamics
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
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

#ifndef PSND_APPL_H
#define PSND_APPL_H

#include "psnd/DataSet.h"
#include "psnd/RuleSet.h"

namespace PROJECT_NS {

class Appl : public std::enable_shared_from_this<Appl> {
    Appl();
    virtual std::shared_ptr<RuleSet> GenerateRuleSet(std::shared_ptr<DataSet> DS);
};
};  // namespace PROJECT_NS

#endif  // PSND_APPL_H

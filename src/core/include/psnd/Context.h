/**@file        Context.h
 * @brief       Declaration of the Context class
 * @details     This file helps in manage context of Kernels
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
 * <tr><td> 2024-04-23  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef PSND_CONTEXT_H
#define PSND_CONTEXT_H

#include "psnd/DataSet.h"
#include "psnd/Model.h"
#include "psnd/Param.h"
#include "psnd/Platform.h"
#include "psnd/RuleSet.h"
#include "psnd/Solver.h"
#include "psnd/System.h"

namespace PROJECT_NS {

class Context {
   public:
    Context(std::shared_ptr<Platform> plat, std::shared_ptr<System> sys, std::vector<std::shared_ptr<Solver>> solvers);

    inline std::shared_ptr<Param>                            getParam() { return _param; }
    inline std::shared_ptr<DataSet>                          getDataSet() { return _dataset; }
    inline std::shared_ptr<RuleSet>                          getRuleSet() { return _ruleset; }
    inline std::shared_ptr<System>                           getSystem() { return _system; }
    inline std::vector<std::vector<std::shared_ptr<Solver>>> getSolvers() { return _solvers; }
    inline std::shared_ptr<Platform>                         getPlatfrom() { return _platform; }

    Status&     run(Status& stat);
    Status&     run_all(Status& stat);
    std::string summary(Status& stat);
    Result      result(Status& stat);

   private:
    std::shared_ptr<Status>                           _stat;
    std::shared_ptr<Param>                            _param;
    std::shared_ptr<DataSet>                          _dataset;
    std::shared_ptr<RuleSet>                          _ruleset;
    std::shared_ptr<System>                           _system;
    std::vector<std::vector<std::shared_ptr<Solver>>> _solvers;
    std::shared_ptr<Platform>                         _platform;
};
};  // namespace PROJECT_NS

#endif  // PSND_CONTEXT_H

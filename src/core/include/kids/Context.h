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
 * <tr><td> 2024-04-23  <td> initial version.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Context_H
#define KIDS_Context_H

#include "kids/DataSet.h"
#include "kids/Model.h"
#include "kids/Param.h"
#include "kids/Platform.h"
#include "kids/RuleSet.h"
#include "kids/Solver.h"
#include "kids/System.h"

namespace PROJECT_NS {

class Context {
    std::shared_ptr<Param>               getParam() { return _param; }
    std::shared_ptr<DataSet>             getDataSet() { return _dataset; }
    std::shared_ptr<RuleSet>             getRuleSet() { return _ruleset; }
    std::shared_ptr<System>              getSystem() { return _system; }
    std::shared_ptr<Model>               getModel() { return _model; }
    std::vector<std::shared_ptr<Solver>> getSolvers() { return _solvers; }
    std::shared_ptr<Platform>            getPlatfrom() { return _platform; }

    Status execute() {
        auto   begin = std::chrono::steady_clock::now();
        Status stat;
        {
            for (auto& solver : _solvers) solver->setInputParam(_param);
            for (auto& solver : _solvers) solver->setInputDataSet(_dataset);
            solver->initializeKernel(stat);  // @necessary?

            // get Monte Carlo Dimension from Param
            int         M         = PM->get_int("M", LOC(), 1);
            std::string directory = PM->get_string("directory", LOC(), "default");
            MPI_Guard   guard(M);
            MPI_Barrier(MPI_COMM_WORLD);
            for (int icalc = guard.istart; icalc < guard.iend; ++icalc) {
                stat.icalc = icalc;
                solver->initializeKernel(stat);
                solver->executeKernel(stat);
                solver->finalizeKernel(stat);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            auto collect = solver->getRuleSet()->getResult(1).data();
            auto reduced = solver->getRuleSet()->getResult(2).data();
            for (int i = 0; i < collect.size(); ++i) {
                std::cout << std::get<0>(collect[i]) << "\n";
                std::cout << std::get<0>(reduced[i]) << "\n";
                auto [key1, from_data, type1, size1, nframe1] = collect[i];
                auto [key2, to_data, type2, size2, nframe2]   = reduced[i];
                MPI_Guard::reduce(std::make_tuple(type1, from_data, to_data, size1));
            }
            // report time cost
            if (MPI_Guard::isroot) RuleSet::flush_all(directory, 2);
        }
        auto   end        = std::chrono::steady_clock::now();
        double total_time = static_cast<std::chrono::duration<double>>(end - begin).count();
    }

   private:
    std::shared_ptr<Param>               _param;
    std::shared_ptr<DataSet>             _dataset;
    std::shared_ptr<RuleSet>             _ruleset;
    std::shared_ptr<System>              _system;
    std::shared_ptr<Model>               _model;
    std::vector<std::shared_ptr<Solver>> _solver;
    std::shared_ptr<Platform>            _platform;
};
};  // namespace PROJECT_NS

#endif  // KIDS_Context_H

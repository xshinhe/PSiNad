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

namespace PROJECT_NS {

class Context {
    std::shared_ptr<Param>    getParam() { return _param; }
    std::shared_ptr<DataSet>  getDataSet() { return _dataset; }
    std::shared_ptr<System>   getSystem() { return _system; }
    std::shared_ptr<Model>    getModel() { return _model; }
    std::shared_ptr<Solver>   getSolver() { return _solver; }
    std::shared_ptr<Platform> getPlatfrom() { return _platform; }

    int run();

   private:
    std::shared_ptr<Kernel> find_result() {
        // linear search
        shared_ptr<Kernel> current_ker;
        shared_ptr<Kernel> parent_ker;
        while () {};
    }

    std::shared_ptr<Param>    _param;
    std::shared_ptr<DataSet>  _dataset;
    std::shared_ptr<System>   _system;
    std::shared_ptr<Model>    _model;
    std::shared_ptr<Solver>   _solver;
    std::shared_ptr<Platform> _platform;
};
};  // namespace PROJECT_NS

#endif  // KIDS_Context_H

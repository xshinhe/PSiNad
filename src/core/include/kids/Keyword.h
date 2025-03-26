/**@file        Keyword.h
 * @brief       this file provide mapping Kernel in Param
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of Xin's PhD research conducted by Professor
 *  Liu's Group at the College of Chemistry and Molecular Engineering, Peking
 *  University. All rights are reserved by Peking University. You should have
 *  received a copy of the GNU Lesser General Public License along with this
 *  software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-14  <td> Update the file.
 * </table>
 *
 **********************************************************************************
 */

#ifndef Keyword_H
#define Keyword_H

#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "kids/Param.h"

using DressingParam_functype = std::function<void(std::shared_ptr<Param>)>;

extern std::vector<DressingParam_functype> DressingParam_funcvector;

inline void appl();

void DressingParam_func1(std::shared_ptr<Param> PM) {
    if (!PM->has_key({"solver.scheme"})) return;

    int debug_level = PM->get_int({"debug_level"}, LOC(), 0);

    std::string scheme = PM->get_string({"solver.scheme"}, LOC(), "");
    if (scheme == "") {  // pass
        throw kids_error("blank solver scheme");
    } else if (scheme == "CMM" || scheme == "CPS" || scheme == "CMM/CPS" || scheme == "MF-CMM" ||
               scheme == "MF-CMM/CPS") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Constraint");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.naforce", "EHR");
        PM->set_int_ifndef("solver.hopping_choose_type", 0);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_real_ifndef("solver.gamma", -1.0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "CMMcv" || scheme == "CMMcv/CPS" || scheme == "MFcv-CMM" || scheme == "MFcv-CMM/CPS") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Constraint");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.naforce", "EHR");
        PM->set_int_ifndef("solver.hopping_choose_type", 0);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_real_ifndef("solver.gamma", -1.0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", true);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "Focus" || scheme == "CFS" || scheme == "CMM/CFS" || scheme == "MF-CMM/CFS") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Focus");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.naforce", "EHR");
        PM->set_int_ifndef("solver.hopping_choose_type", 0);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_real_ifndef("solver.gamma", -1.0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "EHR" || scheme == "Ehrenfest" || scheme == "MF-RHO/EHR") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Focus");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.naforce", "EHR");
        PM->set_int_ifndef("solver.hopping_choose_type", 0);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_real_ifndef("solver.gamma", 0.0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "FSSH" || scheme == "SH-RHO/EHR") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Focus");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.nuc_repr_flag", "Adiabatic");
        PM->set_string_ifndef("solver.naforce", "BO");
        PM->set_int_ifndef("solver.hopping_choose_type", 2);
        PM->set_int_ifndef("solver.hopping_direction_type", 1);
        PM->set_real_ifndef("solver.gamma", 0.0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", true);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "MASH") {  // multistate surface hopping
        if (debug_level < 3) throw kids_error("unsupported scheme!");

        PM->set_string_ifndef("solver.sampling_ele_flag", "Constraint");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.nuc_repr_flag", "Adiabatic");
        PM->set_string_ifndef("solver.naforce", "BO");
        PM->set_int_ifndef("solver.hopping_choose_type", 1);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_real_ifndef("solver.gamma", -2.0);
        PM->set_bool_ifndef("solver.reflect", true);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", false);
        PM->set_bool_ifndef("solver.conserve_scale", false);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "NAF" || scheme == "NAF-CMM" || scheme == "NAF-CMM/CPS") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "Constraint");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.nuc_repr_flag", "Adiabatic");
        PM->set_string_ifndef("solver.naforce", "EHR");
        PM->set_int_ifndef("solver.hopping_choose_type", 0);
        PM->set_int_ifndef("solver.hopping_direction_type", 0);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", true);
        PM->set_bool_ifndef("solver.conserve_scale", true);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else if (scheme == "NAFTW" || scheme == "NAF-TW" || scheme == "NAF-TW/TWS") {
        PM->set_string_ifndef("solver.sampling_ele_flag", "SQCtest01");
        PM->set_string_ifndef("solver.sampling_nuc_flag", "Gaussian");
        PM->set_string_ifndef("solver.nuc_repr_flag", "Adiabatic");
        PM->set_string_ifndef("solver.naforce", "NAF");
        PM->set_int_ifndef("solver.hopping_choose_type", 1);
        PM->set_int_ifndef("solver.hopping_direction_type", 2);
        PM->set_bool_ifndef("solver.reflect", false);
        PM->set_bool_ifndef("solver.use_cv", false);
        PM->set_bool_ifndef("solver.use_fssh", false);
        PM->set_bool_ifndef("solver.offd_projected", true);
        PM->set_bool_ifndef("solver.conserve_scale", true);
        PM->set_bool_ifndef("solver.basis_swith", false);
    } else {
        throw kids_error("unsupported scheme!");
    }
    return;
}

// this file

void appl(std::shared_ptr<Param> PM);

#endif  // Keyword_H

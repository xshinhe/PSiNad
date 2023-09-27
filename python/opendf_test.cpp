
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../models/bath/bath.h"
#include "../models/bo_forcefield/liquidne_model.h"
#include "../models/bo_forcefield/md1d_models.h"
#include "../models/bo_forcefield/sctest_models.h"
#include "../models/bo_forcefield/smallmol_models.h"
#include "../models/bo_forcefield/water_models.h"
#include "../models/forcefieldbase.h"
#include "../models/interface/interf_gausstddft.h"
#include "../models/interface/interf_mndo99mrci.h"
#include "../models/model.h"
#include "../models/nad_forcefield/ZnPc.h"
#include "../models/nad_forcefield/atomced.h"
#include "../models/nad_forcefield/lvcm_model.h"
#include "../models/nad_forcefield/manysite_models.h"
#include "../models/nad_forcefield/nad1d_models.h"
#include "../models/nad_forcefield/pyrincavity_models.h"
#include "../models/nad_forcefield/scatter1d_models.h"
#include "../models/nad_forcefield/spectrum_nadmodels.h"
#include "../models/nad_forcefield/systembath.h"
#include "../models/py_interface/interf_pythonff.h"
#include "../models/thermostat/thermostat.h"
#include "../solvers/solver.h"
#include "../solvers/solvers_PI/cmd_solver.h"
#include "../solvers/solvers_PI/pild_solver.h"
#include "../solvers/solvers_el/basis_set.h"
#include "../solvers/solvers_el/simple_integral_gto.h"
#include "../solvers/solvers_el/solver_scf.h"
#include "../solvers/solvers_md/mbpimd_solver.h"
#include "../solvers/solvers_md/mespimd_solver.h"
#include "../solvers/solvers_md/pimd_solver.h"
#include "../solvers/solvers_md/pimdpara_solver.h"
#include "../solvers/solvers_md/ppimd_solver.h"
#include "../solvers/solvers_md/traj.h"
#include "../solvers/solvers_nad/multi_nadtraj.h"
#include "../solvers/solvers_nad/nadtcfer.h"
#include "../solvers/solvers_nad/nadtraj.h"
#include "../solvers/solvers_nad/solver_cmm.h"
#include "../solvers/solvers_nad/solver_lsc.h"
#include "../solvers/solvers_nad/solver_mmd.h"
#include "../solvers/solvers_nad/solver_pmm.h"
#include "../solvers/solvers_nad/solver_prodmps.h"
#include "../solvers/solvers_nad/solver_qcpi.h"
#include "../solvers/solvers_nad/solver_sh.h"
#include "../solvers/solvers_nad/solver_smm.h"
#include "../solvers/solvers_nad/solver_sqc.h"
#include "../solvers/solvers_nad/solver_wmm.h"
#include "../solvers/solvers_rd/solver_heom.h"
#include "../solvers/solvers_rd/solver_redfield.h"
#include "../solvers/solvers_rd/solver_sse.h"
#include "../utils/definitions.h"

namespace py = pybind11;

// clang-format off
PYBIND11_MODULE(libopendf, m) {

    #include "opendf_phys.bind"

    py::module models_m = m.def_submodule("models");
    py::module solvers_m = m.def_submodule("solvers");



    class PyTrampoline_Model : public Model {
        public:
        using Model::Model;
    };

    py::class_<Model, PyTrampoline_Model>(models_m, "Model", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("tag", &Model::tag)
    .def("ref_workr", &Model::ref_workr, py::return_value_policy::reference_internal)
    .def("ref_workc", &Model::ref_workc, py::return_value_policy::reference_internal);

    class PyTrampoline_ForceField : public ForceField {
        public:
        using ForceField::ForceField;
    };

    py::class_<ForceField, Model, PyTrampoline_ForceField>(models_m, "ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("type", &ForceField::type);

    class PyTrampoline_BO_ForceField : public BO_ForceField {
        public:
        using BO_ForceField::BO_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<BO_ForceField, ForceField, PyTrampoline_BO_ForceField>(models_m, "BO_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_mod_M", &BO_ForceField::ref_mod_M, py::return_value_policy::reference_internal)
    .def("ref_mod_W", &BO_ForceField::ref_mod_W, py::return_value_policy::reference_internal)
    .def("ref_mod_R0", &BO_ForceField::ref_mod_R0, py::return_value_policy::reference_internal)
    .def("ref_mod_P0", &BO_ForceField::ref_mod_P0, py::return_value_policy::reference_internal)
    .def("ref_mod_sigmaR", &BO_ForceField::ref_mod_sigmaR, py::return_value_policy::reference_internal)
    .def("ref_mod_sigmaP", &BO_ForceField::ref_mod_sigmaP, py::return_value_policy::reference_internal)
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_Water_ForceField : public Water_ForceField {
        public:
        using Water_ForceField::Water_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Water_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Water_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Water_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<Water_ForceField, BO_ForceField, PyTrampoline_Water_ForceField>(models_m, "Water_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_charge_arr", &Water_ForceField::ref_charge_arr, py::return_value_policy::reference_internal)
    .def("ref_pbox", &Water_ForceField::ref_pbox, py::return_value_policy::reference_internal)
    .def_static("name", &Water_ForceField::name)
    .def("ForceField_init", [](Water_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& icycle) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, icycle); 
        }
    )
    .def("ForceField_spec", [](Water_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("ForceField_npes", [](Water_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_SmallMol_ForceField : public SmallMol_ForceField {
        public:
        using SmallMol_ForceField::SmallMol_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SmallMol_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SmallMol_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SmallMol_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<SmallMol_ForceField, BO_ForceField, PyTrampoline_SmallMol_ForceField>(models_m, "SmallMol_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &SmallMol_ForceField::name)
    .def("ForceField_init", [](SmallMol_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& icycle) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, icycle); 
        }
    )
    .def("ForceField_spec", [](SmallMol_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("ForceField_npes", [](SmallMol_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_SCTEST_ForceField : public SCTEST_ForceField {
        public:
        using SCTEST_ForceField::SCTEST_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SCTEST_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SCTEST_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SCTEST_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<SCTEST_ForceField, BO_ForceField, PyTrampoline_SCTEST_ForceField>(models_m, "SCTEST_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &SCTEST_ForceField::name)
    .def("ForceField_init", [](SCTEST_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& icycle) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, icycle); 
        }
    )
    .def("ForceField_spec", [](SCTEST_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("ForceField_npes", [](SCTEST_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes_SC1D", [](SCTEST_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes_SC1D(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes_SC2D", [](SCTEST_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes_SC2D(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_LiquidNe_ForceField : public LiquidNe_ForceField {
        public:
        using LiquidNe_ForceField::LiquidNe_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LiquidNe_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LiquidNe_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LiquidNe_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<LiquidNe_ForceField, BO_ForceField, PyTrampoline_LiquidNe_ForceField>(models_m, "LiquidNe_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def_static("name", &LiquidNe_ForceField::name)
    .def("ForceField_init", [](LiquidNe_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& icycle) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, icycle); 
        }
    )
    .def("ForceField_spec", [](LiquidNe_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("ForceField_npes", [](LiquidNe_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_MD1D_ForceField : public MD1D_ForceField {
        public:
        using MD1D_ForceField::MD1D_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MD1D_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MD1D_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MD1D_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<MD1D_ForceField, BO_ForceField, PyTrampoline_MD1D_ForceField>(models_m, "MD1D_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &MD1D_ForceField::name)
    .def("ForceField_init", [](MD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& icycle) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, icycle); 
        }
    )
    .def("ForceField_spec", [](MD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("ForceField_npes", [](MD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_Nad_ForceField : public Nad_ForceField {
        public:
        using Nad_ForceField::Nad_ForceField;

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<Nad_ForceField, BO_ForceField, PyTrampoline_Nad_ForceField>(models_m, "Nad_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_mod_eac", &Nad_ForceField::ref_mod_eac, py::return_value_policy::reference_internal)
    .def("ref_mod_rho", &Nad_ForceField::ref_mod_rho, py::return_value_policy::reference_internal)
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_SystemBath_ForceField : public SystemBath_ForceField {
        public:
        using SystemBath_ForceField::SystemBath_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_epes_SpinBoson(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                          const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_SpinBoson, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_epes_SiteExciton(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                            const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_SiteExciton, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_epes_General(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_General, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int get_nbath() override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            get_nbath, // func name
            
            );
        }

        int get_Nb() override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            get_Nb, // func name
            
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<SystemBath_ForceField, Nad_ForceField, PyTrampoline_SystemBath_ForceField>(models_m, "SystemBath_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("L", &SystemBath_ForceField::L)
    .def("ref_Hsys", &SystemBath_ForceField::ref_Hsys, py::return_value_policy::reference_internal)
    .def("ref_Q", &SystemBath_ForceField::ref_Q, py::return_value_policy::reference_internal)
    .def("ref_omegas", &SystemBath_ForceField::ref_omegas, py::return_value_policy::reference_internal)
    .def("ref_coeffs", &SystemBath_ForceField::ref_coeffs, py::return_value_policy::reference_internal)
    .def("ref_CL", &SystemBath_ForceField::ref_CL, py::return_value_policy::reference_internal)
    .def("ref_QL", &SystemBath_ForceField::ref_QL, py::return_value_policy::reference_internal)
    .def("ref_Xnj", &SystemBath_ForceField::ref_Xnj, py::return_value_policy::reference_internal)
    .def_static("name", &SystemBath_ForceField::name)
    .def("ForceField_npes", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_SpinBoson", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_SpinBoson(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_SiteExciton", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_SiteExciton(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_General", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_General(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("get_nbath", static_cast<int (SystemBath_ForceField::*)()>(&SystemBath_ForceField::get_nbath))
    .def("get_Nb", static_cast<int (SystemBath_ForceField::*)()>(&SystemBath_ForceField::get_Nb))
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_NAD1D_ForceField : public NAD1D_ForceField {
        public:
        using NAD1D_ForceField::NAD1D_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int NAD1D_plot() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            NAD1D_plot, // func name
            
            );
        }

        int ForceField_epes_Morse3A(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_Morse3A, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_Morse3B(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_Morse3B, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_Morse3C(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_Morse3C, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_Morse15(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_Morse15, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_IVP1(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_IVP1, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_IVP2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_IVP2, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_IVP3(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_IVP3, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_IVP4(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_IVP4, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_CL1D(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_CL1D, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_JC1D(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_JC1D, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_NA_I(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NAD1D_ForceField, // parent class
            ForceField_epes_NA_I, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<NAD1D_ForceField, Nad_ForceField, PyTrampoline_NAD1D_ForceField>(models_m, "NAD1D_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &NAD1D_ForceField::name)
    .def("ForceField_spec", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("ForceField_npes", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("NAD1D_plot", static_cast<int (NAD1D_ForceField::*)()>(&NAD1D_ForceField::NAD1D_plot))
    .def("ForceField_epes_Morse3A", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_Morse3A(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_Morse3B", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_Morse3B(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_Morse3C", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_Morse3C(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_Morse15", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_Morse15(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_IVP1", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_IVP1(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_IVP2", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_IVP2(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_IVP3", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_IVP3(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_IVP4", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_IVP4(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_CL1D", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_CL1D(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_JC1D", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_JC1D(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_NA_I", [](NAD1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_NA_I(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_Spectrum_NAD_ForceField : public Spectrum_NAD_ForceField {
        public:
        using Spectrum_NAD_ForceField::Spectrum_NAD_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Spectrum_NAD_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Spectrum_NAD_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Spectrum_NAD_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<Spectrum_NAD_ForceField, Nad_ForceField, PyTrampoline_Spectrum_NAD_ForceField>(models_m, "Spectrum_NAD_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("first_call", &Spectrum_NAD_ForceField::first_call)
    .def_readwrite("Fminus1", &Spectrum_NAD_ForceField::Fminus1)
    .def_readwrite("ground_shift", &Spectrum_NAD_ForceField::ground_shift)
    .def("ref_workr_v", &Spectrum_NAD_ForceField::ref_workr_v, py::return_value_policy::reference_internal)
    .def("ref_workr_dv", &Spectrum_NAD_ForceField::ref_workr_dv, py::return_value_policy::reference_internal)
    .def_static("name", &Spectrum_NAD_ForceField::name)
    .def("ForceField_npes", [](Spectrum_NAD_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](Spectrum_NAD_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_LVCM_ForceField : public LVCM_ForceField {
        public:
        using LVCM_ForceField::LVCM_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_PYR(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_epes_PYR, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_CrCO5_2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_epes_CrCO5_2, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_CrCO5_5(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_epes_CrCO5_5, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_BEN_5(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                      const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LVCM_ForceField, // parent class
            ForceField_epes_BEN_5, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<LVCM_ForceField, Nad_ForceField, PyTrampoline_LVCM_ForceField>(models_m, "LVCM_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_Hsys", &LVCM_ForceField::ref_Hsys, py::return_value_policy::reference_internal)
    .def("ref_ECI", &LVCM_ForceField::ref_ECI, py::return_value_policy::reference_internal)
    .def("ref_KCI", &LVCM_ForceField::ref_KCI, py::return_value_policy::reference_internal)
    .def_static("name", &LVCM_ForceField::name)
    .def("ForceField_npes", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_PYR", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_PYR(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_CrCO5_2", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_CrCO5_2(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_CrCO5_5", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_CrCO5_5(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_BEN_5", [](LVCM_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_BEN_5(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_ManySite_ForceField : public ManySite_ForceField {
        public:
        using ManySite_ForceField::ManySite_ForceField;

        int ForceField_heff(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ManySite_ForceField, // parent class
            ForceField_heff, // func name
            H, rhos, mdim, fdim
            );
        }

        int ForceField_heff_Ising(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ManySite_ForceField, // parent class
            ForceField_heff_Ising, // func name
            H, rhos, mdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<ManySite_ForceField, Nad_ForceField, PyTrampoline_ManySite_ForceField>(models_m, "ManySite_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("M", &ManySite_ForceField::M)
    .def_readwrite("Jp", &ManySite_ForceField::Jp)
    .def_readwrite("Jz", &ManySite_ForceField::Jz)
    .def_readwrite("alpha", &ManySite_ForceField::alpha)
    .def_readwrite("omega", &ManySite_ForceField::omega)
    .def("ref_JpMat", &ManySite_ForceField::ref_JpMat, py::return_value_policy::reference_internal)
    .def("ref_JzMat", &ManySite_ForceField::ref_JzMat, py::return_value_policy::reference_internal)
    .def("ref_redX", &ManySite_ForceField::ref_redX, py::return_value_policy::reference_internal)
    .def("ref_redY", &ManySite_ForceField::ref_redY, py::return_value_policy::reference_internal)
    .def("ref_redZ", &ManySite_ForceField::ref_redZ, py::return_value_policy::reference_internal)
    .def_static("name", &ManySite_ForceField::name)
    .def("get_M", static_cast<int (ManySite_ForceField::*)()>(&ManySite_ForceField::get_M))
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_Scatter1D_ForceField : public Scatter1D_ForceField {
        public:
        using Scatter1D_ForceField::Scatter1D_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int Scatter1D_plot(const double& Xrange) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            Scatter1D_plot, // func name
            Xrange
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_SAC(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_SAC, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_SAC2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_SAC2, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_DAC(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_DAC, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_ECR(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_ECR, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_DBG(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_DBG, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_DAG(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_DAG, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_DRN(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Scatter1D_ForceField, // parent class
            ForceField_epes_DRN, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<Scatter1D_ForceField, Nad_ForceField, PyTrampoline_Scatter1D_ForceField>(models_m, "Scatter1D_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &Scatter1D_ForceField::name)
    .def("ForceField_spec", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("ForceField_npes", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("Scatter1D_plot", static_cast<int (Scatter1D_ForceField::*)(const double&)>(&Scatter1D_ForceField::Scatter1D_plot))
    .def("ForceField_epes", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_SAC", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_SAC(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_SAC2", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_SAC2(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_DAC", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_DAC(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_ECR", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_ECR(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_DBG", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_DBG(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_DAG", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_DAG(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_DRN", [](Scatter1D_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_DRN(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_PyrCav_ForceField : public PyrCav_ForceField {
        public:
        using PyrCav_ForceField::PyrCav_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_PC1(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_epes_PC1, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_PC2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_epes_PC2, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_epes_PC3(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PyrCav_ForceField, // parent class
            ForceField_epes_PC3, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<PyrCav_ForceField, Nad_ForceField, PyTrampoline_PyrCav_ForceField>(models_m, "PyrCav_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_WCI", &PyrCav_ForceField::ref_WCI, py::return_value_policy::reference_internal)
    .def("ref_ECI", &PyrCav_ForceField::ref_ECI, py::return_value_policy::reference_internal)
    .def("ref_KCI", &PyrCav_ForceField::ref_KCI, py::return_value_policy::reference_internal)
    .def_static("name", &PyrCav_ForceField::name)
    .def("ForceField_npes", [](PyrCav_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](PyrCav_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_PC1", [](PyrCav_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_PC1(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_PC2", [](PyrCav_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_PC2(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes_PC3", [](PyrCav_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes_PC3(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_AtomCED_ForceField : public AtomCED_ForceField {
        public:
        using AtomCED_ForceField::AtomCED_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            AtomCED_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            AtomCED_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            AtomCED_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            AtomCED_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int get_Nb() override {
            PYBIND11_OVERRIDE(
            int, // return type
            AtomCED_ForceField, // parent class
            get_Nb, // func name
            
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<AtomCED_ForceField, Nad_ForceField, PyTrampoline_AtomCED_ForceField>(models_m, "AtomCED_ForceField", py::dynamic_attr())
    .def(py::init<const Param&, const int&>())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_Hsys", &AtomCED_ForceField::ref_Hsys, py::return_value_policy::reference_internal)
    .def("ref_omegas", &AtomCED_ForceField::ref_omegas, py::return_value_policy::reference_internal)
    .def("ref_coeffs", &AtomCED_ForceField::ref_coeffs, py::return_value_policy::reference_internal)
    .def_static("name", &AtomCED_ForceField::name)
    .def("ForceField_npes", [](AtomCED_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](AtomCED_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_Nb", static_cast<int (AtomCED_ForceField::*)()>(&AtomCED_ForceField::get_Nb))
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_ZnPc_ForceField : public ZnPc_ForceField {
        public:
        using ZnPc_ForceField::ZnPc_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ZnPc_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ZnPc_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ZnPc_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int get_nbath() override {
            PYBIND11_OVERRIDE(
            int, // return type
            ZnPc_ForceField, // parent class
            get_nbath, // func name
            
            );
        }

        int get_Nb() override {
            PYBIND11_OVERRIDE(
            int, // return type
            ZnPc_ForceField, // parent class
            get_Nb, // func name
            
            );
        }

        int ForceField_epes_SpinBoson(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                          const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_SpinBoson, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_epes_SiteExciton(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                            const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_SiteExciton, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_epes_General(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            ForceField_epes_General, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SystemBath_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<ZnPc_ForceField, SystemBath_ForceField, PyTrampoline_ZnPc_ForceField>(models_m, "ZnPc_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_idxarr_L", &ZnPc_ForceField::ref_idxarr_L, py::return_value_policy::reference_internal)
    .def("ref_idxarr_es", &ZnPc_ForceField::ref_idxarr_es, py::return_value_policy::reference_internal)
    .def("ref_idxarr_hs", &ZnPc_ForceField::ref_idxarr_hs, py::return_value_policy::reference_internal)
    .def("ref_dQe1", &ZnPc_ForceField::ref_dQe1, py::return_value_policy::reference_internal)
    .def("ref_dQe2", &ZnPc_ForceField::ref_dQe2, py::return_value_policy::reference_internal)
    .def("ref_dQc", &ZnPc_ForceField::ref_dQc, py::return_value_policy::reference_internal)
    .def("ref_dQa", &ZnPc_ForceField::ref_dQa, py::return_value_policy::reference_internal)
    .def("ref_w2dQe1", &ZnPc_ForceField::ref_w2dQe1, py::return_value_policy::reference_internal)
    .def("ref_w2dQe2", &ZnPc_ForceField::ref_w2dQe2, py::return_value_policy::reference_internal)
    .def("ref_w2dQc", &ZnPc_ForceField::ref_w2dQc, py::return_value_policy::reference_internal)
    .def("ref_w2dQa", &ZnPc_ForceField::ref_w2dQa, py::return_value_policy::reference_internal)
    .def("ref_Etilde", &ZnPc_ForceField::ref_Etilde, py::return_value_policy::reference_internal)
    .def("ref_Vtilde", &ZnPc_ForceField::ref_Vtilde, py::return_value_policy::reference_internal)
    .def("ref_te_tilde", &ZnPc_ForceField::ref_te_tilde, py::return_value_policy::reference_internal)
    .def("ref_th_tilde", &ZnPc_ForceField::ref_th_tilde, py::return_value_policy::reference_internal)
    .def("ref_tect_tilde", &ZnPc_ForceField::ref_tect_tilde, py::return_value_policy::reference_internal)
    .def("ref_thct_tilde", &ZnPc_ForceField::ref_thct_tilde, py::return_value_policy::reference_internal)
    .def("ref_eigen_E", &ZnPc_ForceField::ref_eigen_E, py::return_value_policy::reference_internal)
    .def("ref_eigen_T", &ZnPc_ForceField::ref_eigen_T, py::return_value_policy::reference_internal)
    .def_static("name", &ZnPc_ForceField::name)
    .def("index", static_cast<int (ZnPc_ForceField::*)(const int&, const int&, const int&)>(&ZnPc_ForceField::index))
    .def("init_Hamiltonian", static_cast<int (ZnPc_ForceField::*)()>(&ZnPc_ForceField::init_Hamiltonian))
    .def("ForceField_npes", [](ZnPc_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](ZnPc_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("get_nbath", static_cast<int (ZnPc_ForceField::*)()>(&ZnPc_ForceField::get_nbath))
    .def("get_Nb", static_cast<int (ZnPc_ForceField::*)()>(&ZnPc_ForceField::get_Nb))
    .def("ForceField_npes", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_SpinBoson", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_SpinBoson(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_SiteExciton", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_SiteExciton(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes_General", [](SystemBath_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes_General(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("get_nbath", static_cast<int (SystemBath_ForceField::*)()>(&SystemBath_ForceField::get_nbath))
    .def("get_Nb", static_cast<int (SystemBath_ForceField::*)()>(&SystemBath_ForceField::get_Nb))
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_GAUSS16_ForceField : public GAUSS16_ForceField {
        public:
        using GAUSS16_ForceField::GAUSS16_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle) override {
            PYBIND11_OVERRIDE(
            int, // return type
            GAUSS16_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle
            );
        }

        int ForceField_epes(num_real* E, num_real* dE, num_real* ddE,
                                num_real* R,  
                                const int& flag, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            GAUSS16_ForceField, // parent class
            ForceField_epes, // func name
            E, dE, ddE, R, flag, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<GAUSS16_ForceField, Nad_ForceField, PyTrampoline_GAUSS16_ForceField>(models_m, "GAUSS16_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_atoms", &GAUSS16_ForceField::ref_atoms, py::return_value_policy::reference_internal)
    .def("ref_mod_Hess", &GAUSS16_ForceField::ref_mod_Hess, py::return_value_policy::reference_internal)
    .def("ref_mod_Tmat", &GAUSS16_ForceField::ref_mod_Tmat, py::return_value_policy::reference_internal)
    .def("ref_nr_samp", &GAUSS16_ForceField::ref_nr_samp, py::return_value_policy::reference_internal)
    .def("ref_np_samp", &GAUSS16_ForceField::ref_np_samp, py::return_value_policy::reference_internal)
    .def_static("name", &GAUSS16_ForceField::name)
    .def("ForceField_epes", [](GAUSS16_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> E_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dE_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddE_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(E_arr.mutable_data(), dE_arr.mutable_data(), ddE_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("parse_g16", static_cast<int (GAUSS16_ForceField::*)(const std::string&)>(&GAUSS16_ForceField::parse_g16))
    .def("calc_hess", [](GAUSS16_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& rdim) {
            return self.calc_hess(R_arr.mutable_data(), rdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_MNDO99_ForceField : public MNDO99_ForceField {
        public:
        using MNDO99_ForceField::MNDO99_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MNDO99_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj
            );
        }

        int ForceField_epes(num_real* E, num_real* dE, num_real* ddE,
                                num_real* R,                       
                                const int& flag,                   
                                const int& rdim, const int& fdim,  
                                const int& itraj, const int& istep) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MNDO99_ForceField, // parent class
            ForceField_epes, // func name
            E, dE, ddE, R, flag, rdim, fdim, itraj, istep
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& istep) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MNDO99_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, istep
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<MNDO99_ForceField, Nad_ForceField, PyTrampoline_MNDO99_ForceField>(models_m, "MNDO99_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def("ref_atoms", &MNDO99_ForceField::ref_atoms, py::return_value_policy::reference_internal)
    .def("ref_mod_Hess", &MNDO99_ForceField::ref_mod_Hess, py::return_value_policy::reference_internal)
    .def("ref_mod_Tmat", &MNDO99_ForceField::ref_mod_Tmat, py::return_value_policy::reference_internal)
    .def("ref_nr_samp", &MNDO99_ForceField::ref_nr_samp, py::return_value_policy::reference_internal)
    .def("ref_np_samp", &MNDO99_ForceField::ref_np_samp, py::return_value_policy::reference_internal)
    .def("ref_nac_prev", &MNDO99_ForceField::ref_nac_prev, py::return_value_policy::reference_internal)
    .def("ref_nac", &MNDO99_ForceField::ref_nac, py::return_value_policy::reference_internal)
    .def("ref_ener", &MNDO99_ForceField::ref_ener, py::return_value_policy::reference_internal)
    .def("ref_ener2", &MNDO99_ForceField::ref_ener2, py::return_value_policy::reference_internal)
    .def("ref_grad", &MNDO99_ForceField::ref_grad, py::return_value_policy::reference_internal)
    .def_static("name", &MNDO99_ForceField::name)
    .def("ForceField_epes", [](MNDO99_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> E_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dE_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddE_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& istep) {
            return self.ForceField_epes(E_arr.mutable_data(), dE_arr.mutable_data(), ddE_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, istep); 
        }
    )
    .def("parse_mndo99", static_cast<int (MNDO99_ForceField::*)(const std::string&)>(&MNDO99_ForceField::parse_mndo99))
    .def("new_keyword", static_cast<std::string (MNDO99_ForceField::*)(const MNDO99KW_map&)>(&MNDO99_ForceField::new_keyword))
    .def("new_task", [](MNDO99_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const std::string& file, 
        const std::string& task_flag, 
        const int& rdim) {
            return self.new_task(R_arr.mutable_data(), file, task_flag, rdim); 
        }
    )
    .def("track_nac_sign", static_cast<int (MNDO99_ForceField::*)()>(&MNDO99_ForceField::track_nac_sign))
    .def("parse_standard", static_cast<int (MNDO99_ForceField::*)(const std::string&)>(&MNDO99_ForceField::parse_standard))
    .def("parse_hessian", static_cast<int (MNDO99_ForceField::*)(const std::string&)>(&MNDO99_ForceField::parse_hessian))
    .def("calc_normalmode", [](MNDO99_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& rdim) {
            return self.calc_normalmode(R_arr.mutable_data(), rdim); 
        }
    )
    .def("calc_samp", static_cast<int (MNDO99_ForceField::*)()>(&MNDO99_ForceField::calc_samp))
    .def("calc_scan", static_cast<int (MNDO99_ForceField::*)()>(&MNDO99_ForceField::calc_scan))
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_PythonFF_ForceField : public PythonFF_ForceField {
        public:
        using PythonFF_ForceField::PythonFF_ForceField;

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PythonFF_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PythonFF_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PythonFF_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim, fdim
            );
        }

        int nspec() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            nspec, // func name
            
            );
        }

        int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_npes, // func name
            V, dV, ddV, R, P, flag, rdim, itraj, isamp
            );
        }

        int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_epes, // func name
            V, dV, ddV, R, flag, rdim, fdim, itraj, isamp
            );
        }

        int CheckForceField() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            CheckForceField, // func name
            
            );
        }

        int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            ForceField_write, // func name
            ofs0, ofs1, nr, np, nm, erho, eeac, eocc, rdim, fdim, itraj, isamp
            );
        }

        int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Nad_ForceField, // parent class
            reduce_force, // func name
            fx, rho, dH, rdim, fdim
            );
        }

        int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_init, // func name
            nr, np, nm, rdim, itraj
            );
        }

        int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_spec, // func name
            nr, np, nm, rdim
            );
        }

        int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp) override {
            PYBIND11_OVERRIDE(
            int, // return type
            BO_ForceField, // parent class
            ForceField_write, // func name
            ofs0, nr, np, nm, rdim, itraj, isamp
            );
        }
    };

    py::class_<PythonFF_ForceField, Nad_ForceField, PyTrampoline_PythonFF_ForceField>(models_m, "PythonFF_ForceField", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_static("name", &PythonFF_ForceField::name)
    .def("ForceField_npes", [](PythonFF_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_epes", [](PythonFF_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("get_F", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::get_F))
    .def("ForceField_spec", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, fdim); 
        }
    )
    .def("nspec", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::nspec))
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim); 
        }
    )
    .def("ForceField_epes", [](Nad_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        const int& flag, 
        const int& rdim, 
        const int& fdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_epes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), flag, rdim, fdim, itraj, isamp); 
        }
    )
    .def("CheckForceField", static_cast<int (Nad_ForceField::*)()>(&Nad_ForceField::CheckForceField))
    .def("get_N", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_N))
    .def("get_Ndim", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::get_Ndim))
    .def("Suggest_dt", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_dt))
    .def("Suggest_tend", static_cast<double (BO_ForceField::*)()>(&BO_ForceField::Suggest_tend))
    .def("ForceField_init", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("ForceField_spec", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim) {
            return self.ForceField_spec(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim); 
        }
    )
    .def("nspec", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::nspec))
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim); 
        }
    )
    .def("ForceField_npes", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> V_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> ddV_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> R_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> P_arr, 
        const int& flag, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_npes(V_arr.mutable_data(), dV_arr.mutable_data(), ddV_arr.mutable_data(), R_arr.mutable_data(), P_arr.mutable_data(), flag, rdim, itraj, isamp); 
        }
    )
    .def("ForceField_init_default_build", static_cast<int (BO_ForceField::*)(const double&, const int&)>(&BO_ForceField::ForceField_init_default_build))
    .def("ForceField_init_default", [](BO_ForceField& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj) {
            return self.ForceField_init_default(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj); 
        }
    )
    .def("CheckForceField", static_cast<int (BO_ForceField::*)()>(&BO_ForceField::CheckForceField))
    .def("ForceField_write", [](BO_ForceField& self, std::ofstream& ofs0, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const int& rdim, 
        const int& itraj, 
        const int& isamp) {
            return self.ForceField_write(ofs0, nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), rdim, itraj, isamp); 
        }
    );

    class PyTrampoline_Solver : public Solver {
        public:
        using Solver::Solver;

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<Solver, Model, PyTrampoline_Solver>(solvers_m, "Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_readwrite("save", &Solver::save)
    .def_readwrite("para_flag", &Solver::para_flag)
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_Traj_Solver : public Traj_Solver {
        public:
        using Traj_Solver::Traj_Solver;

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<Traj_Solver, Solver, PyTrampoline_Traj_Solver>(solvers_m, "Traj_Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_Htot", &Traj_Solver::ref_Htot, py::return_value_policy::reference_internal)
    .def("ref_Ltot", &Traj_Solver::ref_Ltot, py::return_value_policy::reference_internal)
    .def("ref_Stot", &Traj_Solver::ref_Stot, py::return_value_policy::reference_internal)
    .def("ref_Ktot", &Traj_Solver::ref_Ktot, py::return_value_policy::reference_internal)
    .def("ref_Vtot", &Traj_Solver::ref_Vtot, py::return_value_policy::reference_internal)
    .def("ref_nr0", &Traj_Solver::ref_nr0, py::return_value_policy::reference_internal)
    .def("ref_np0", &Traj_Solver::ref_np0, py::return_value_policy::reference_internal)
    .def("ref_nr", &Traj_Solver::ref_nr, py::return_value_policy::reference_internal)
    .def("ref_np", &Traj_Solver::ref_np, py::return_value_policy::reference_internal)
    .def("ref_nm", &Traj_Solver::ref_nm, py::return_value_policy::reference_internal)
    .def("ref_nf", &Traj_Solver::ref_nf, py::return_value_policy::reference_internal)
    .def("ref_vpes", &Traj_Solver::ref_vpes, py::return_value_policy::reference_internal)
    .def("ref_grad", &Traj_Solver::ref_grad, py::return_value_policy::reference_internal)
    .def("ref_hess", &Traj_Solver::ref_hess, py::return_value_policy::reference_internal)
    .def_static("name", &Traj_Solver::name)
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_PIMDTraj_Solver : public PIMDTraj_Solver {
        public:
        using PIMDTraj_Solver::PIMDTraj_Solver;

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_p_harm(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<PIMDTraj_Solver, Traj_Solver, PyTrampoline_PIMDTraj_Solver>(solvers_m, "PIMDTraj_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_nrs", &PIMDTraj_Solver::ref_nrs, py::return_value_policy::reference_internal)
    .def("ref_nks", &PIMDTraj_Solver::ref_nks, py::return_value_policy::reference_internal)
    .def("ref_nps", &PIMDTraj_Solver::ref_nps, py::return_value_policy::reference_internal)
    .def("ref_nms", &PIMDTraj_Solver::ref_nms, py::return_value_policy::reference_internal)
    .def("ref_nfs", &PIMDTraj_Solver::ref_nfs, py::return_value_policy::reference_internal)
    .def("ref_nfks", &PIMDTraj_Solver::ref_nfks, py::return_value_policy::reference_internal)
    .def("ref_masswgt", &PIMDTraj_Solver::ref_masswgt, py::return_value_policy::reference_internal)
    .def("ref_bfwgt", &PIMDTraj_Solver::ref_bfwgt, py::return_value_policy::reference_internal)
    .def("ref_D2", &PIMDTraj_Solver::ref_D2, py::return_value_policy::reference_internal)
    .def("ref_Tran", &PIMDTraj_Solver::ref_Tran, py::return_value_policy::reference_internal)
    .def("ref_vpeses", &PIMDTraj_Solver::ref_vpeses, py::return_value_policy::reference_internal)
    .def("ref_grads", &PIMDTraj_Solver::ref_grads, py::return_value_policy::reference_internal)
    .def("ref_hesses", &PIMDTraj_Solver::ref_hesses, py::return_value_policy::reference_internal)
    .def_static("name", &PIMDTraj_Solver::name)
    .def("ff_calc1", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::init))
    .def("all_X2K", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDTraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("pseudo_inv", [](PIMDTraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_read))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_MESPIMDTraj_Solver : public MESPIMDTraj_Solver {
        public:
        using MESPIMDTraj_Solver::MESPIMDTraj_Solver;

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MESPIMDTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int ff_calc1(const int& level = 1) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MESPIMDTraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int ff_Oi(int i, num_real* Oi_pos, num_real* dVi_pos, num_real* Oi, num_real* dOi, num_real* Vi,
                      num_real* dVi) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MESPIMDTraj_Solver, // parent class
            ff_Oi, // func name
            i, Oi_pos, dVi_pos, Oi, dOi, Vi, dVi
            );
        }

        int ff_OO() override {
            PYBIND11_OVERRIDE(
            int, // return type
            MESPIMDTraj_Solver, // parent class
            ff_OO, // func name
            
            );
        }

        num_real esti_term1(num_real* Q, const bool& dressed = true, const bool& fixed = false) override {
            PYBIND11_OVERRIDE(
            num_real, // return type
            MESPIMDTraj_Solver, // parent class
            esti_term1, // func name
            Q, dressed, fixed
            );
        }

        num_real esti_term2(num_real* Q1, num_real* Q2, const bool& dressed1 = true, const bool& dressed2 = true,
                                const bool& fixed1 = false, const bool& fixed2 = false) override {
            PYBIND11_OVERRIDE(
            num_real, // return type
            MESPIMDTraj_Solver, // parent class
            esti_term2, // func name
            Q1, Q2, dressed1, dressed2, fixed1, fixed2
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MESPIMDTraj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_p_harm(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<MESPIMDTraj_Solver, PIMDTraj_Solver, PyTrampoline_MESPIMDTraj_Solver>(solvers_m, "MESPIMDTraj_Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_Vs", &MESPIMDTraj_Solver::ref_Vs, py::return_value_policy::reference_internal)
    .def("ref_dVs", &MESPIMDTraj_Solver::ref_dVs, py::return_value_policy::reference_internal)
    .def("ref_ddVs", &MESPIMDTraj_Solver::ref_ddVs, py::return_value_policy::reference_internal)
    .def("ref_Es", &MESPIMDTraj_Solver::ref_Es, py::return_value_policy::reference_internal)
    .def("ref_dEs", &MESPIMDTraj_Solver::ref_dEs, py::return_value_policy::reference_internal)
    .def("ref_Ts", &MESPIMDTraj_Solver::ref_Ts, py::return_value_policy::reference_internal)
    .def("ref_mat_fA", &MESPIMDTraj_Solver::ref_mat_fA, py::return_value_policy::reference_internal)
    .def("ref_mat_fD", &MESPIMDTraj_Solver::ref_mat_fD, py::return_value_policy::reference_internal)
    .def("ref_O_pos", &MESPIMDTraj_Solver::ref_O_pos, py::return_value_policy::reference_internal)
    .def("ref_dV_pos", &MESPIMDTraj_Solver::ref_dV_pos, py::return_value_policy::reference_internal)
    .def("ref_O", &MESPIMDTraj_Solver::ref_O, py::return_value_policy::reference_internal)
    .def("ref_OO", &MESPIMDTraj_Solver::ref_OO, py::return_value_policy::reference_internal)
    .def("ref_OOb", &MESPIMDTraj_Solver::ref_OOb, py::return_value_policy::reference_internal)
    .def("ref_OOe", &MESPIMDTraj_Solver::ref_OOe, py::return_value_policy::reference_internal)
    .def("ref_OObe", &MESPIMDTraj_Solver::ref_OObe, py::return_value_policy::reference_internal)
    .def("ref_dO", &MESPIMDTraj_Solver::ref_dO, py::return_value_policy::reference_internal)
    .def("ref_V2", &MESPIMDTraj_Solver::ref_V2, py::return_value_policy::reference_internal)
    .def("ref_hRdOO", &MESPIMDTraj_Solver::ref_hRdOO, py::return_value_policy::reference_internal)
    .def("ref_hRcdOO", &MESPIMDTraj_Solver::ref_hRcdOO, py::return_value_policy::reference_internal)
    .def("ref_hRdOVO", &MESPIMDTraj_Solver::ref_hRdOVO, py::return_value_policy::reference_internal)
    .def("ref_hRcdOVO", &MESPIMDTraj_Solver::ref_hRcdOVO, py::return_value_policy::reference_internal)
    .def("ref_rho", &MESPIMDTraj_Solver::ref_rho, py::return_value_policy::reference_internal)
    .def("ref_rho_op", &MESPIMDTraj_Solver::ref_rho_op, py::return_value_policy::reference_internal)
    .def("ref_eac0", &MESPIMDTraj_Solver::ref_eac0, py::return_value_policy::reference_internal)
    .def("ref_rho0", &MESPIMDTraj_Solver::ref_rho0, py::return_value_policy::reference_internal)
    .def_static("name", &MESPIMDTraj_Solver::name)
    .def("init", static_cast<int (MESPIMDTraj_Solver::*)(const int&)>(&MESPIMDTraj_Solver::init))
    .def("ff_calc1", static_cast<int (MESPIMDTraj_Solver::*)(const int&)>(&MESPIMDTraj_Solver::ff_calc1))
    .def("ff_Oi", [](MESPIMDTraj_Solver& self, int i, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> Oi_pos_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dVi_pos_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> Oi_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dOi_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> Vi_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> dVi_arr) {
            return self.ff_Oi(i, Oi_pos_arr.mutable_data(), dVi_pos_arr.mutable_data(), Oi_arr.mutable_data(), dOi_arr.mutable_data(), Vi_arr.mutable_data(), dVi_arr.mutable_data()); 
        }
    )
    .def("ff_OO", static_cast<int (MESPIMDTraj_Solver::*)()>(&MESPIMDTraj_Solver::ff_OO))
    .def("esti_term1", [](MESPIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> Q_arr, 
        const bool& dressed, 
        const bool& fixed) {
            return self.esti_term1(Q_arr.mutable_data(), dressed, fixed); 
        }
    )
    .def("esti_term2", [](MESPIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> Q1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> Q2_arr, 
        const bool& dressed1, 
        const bool& dressed2, 
        const bool& fixed1, 
        const bool& fixed2) {
            return self.esti_term2(Q1_arr.mutable_data(), Q2_arr.mutable_data(), dressed1, dressed2, fixed1, fixed2); 
        }
    )
    .def("estimator", static_cast<int (MESPIMDTraj_Solver::*)(const int&, TCFnucl&)>(&MESPIMDTraj_Solver::estimator))
    .def("ff_calc1", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::init))
    .def("all_X2K", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDTraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("pseudo_inv", [](PIMDTraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_read))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_PIMDPARATraj_Solver : public PIMDPARATraj_Solver {
        public:
        using PIMDPARATraj_Solver::PIMDPARATraj_Solver;

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_p_harm(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int all_X2K() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            all_X2K, // func name
            
            );
        }

        int all_K2X() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            all_K2X, // func name
            
            );
        }

        int all_FX2FK() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            all_FX2FK, // func name
            
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int cons_rot() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            cons_rot, // func name
            
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDPARATraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<PIMDPARATraj_Solver, Traj_Solver, PyTrampoline_PIMDPARATraj_Solver>(solvers_m, "PIMDPARATraj_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def("ref_nrs", &PIMDPARATraj_Solver::ref_nrs, py::return_value_policy::reference_internal)
    .def("ref_nks", &PIMDPARATraj_Solver::ref_nks, py::return_value_policy::reference_internal)
    .def("ref_nps", &PIMDPARATraj_Solver::ref_nps, py::return_value_policy::reference_internal)
    .def("ref_nms", &PIMDPARATraj_Solver::ref_nms, py::return_value_policy::reference_internal)
    .def("ref_nfs", &PIMDPARATraj_Solver::ref_nfs, py::return_value_policy::reference_internal)
    .def("ref_nfks", &PIMDPARATraj_Solver::ref_nfks, py::return_value_policy::reference_internal)
    .def("ref_masswgt", &PIMDPARATraj_Solver::ref_masswgt, py::return_value_policy::reference_internal)
    .def("ref_bfwgt", &PIMDPARATraj_Solver::ref_bfwgt, py::return_value_policy::reference_internal)
    .def("ref_D2", &PIMDPARATraj_Solver::ref_D2, py::return_value_policy::reference_internal)
    .def("ref_Tran", &PIMDPARATraj_Solver::ref_Tran, py::return_value_policy::reference_internal)
    .def("ref_vpeses", &PIMDPARATraj_Solver::ref_vpeses, py::return_value_policy::reference_internal)
    .def("ref_grads", &PIMDPARATraj_Solver::ref_grads, py::return_value_policy::reference_internal)
    .def("ref_hesses", &PIMDPARATraj_Solver::ref_hesses, py::return_value_policy::reference_internal)
    .def_static("name", &PIMDPARATraj_Solver::name)
    .def("ff_calc1", static_cast<int (PIMDPARATraj_Solver::*)(const int&)>(&PIMDPARATraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDPARATraj_Solver::*)(int&, int)>(&PIMDPARATraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDPARATraj_Solver::*)(int&, int)>(&PIMDPARATraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDPARATraj_Solver::*)(TCFnucl&, const int&)>(&PIMDPARATraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDPARATraj_Solver::*)(TCFnucl&, const int&)>(&PIMDPARATraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDPARATraj_Solver::*)(const num_real&)>(&PIMDPARATraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDPARATraj_Solver::*)(const int&, TCFnucl&)>(&PIMDPARATraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDPARATraj_Solver::*)(const int&, TCFnucl&)>(&PIMDPARATraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDPARATraj_Solver::*)(const int&)>(&PIMDPARATraj_Solver::init))
    .def("final", static_cast<int (PIMDPARATraj_Solver::*)(const int&)>(&PIMDPARATraj_Solver::final))
    .def("all_X2K", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDPARATraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("cons_rot", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::cons_rot))
    .def("pseudo_inv", [](PIMDPARATraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDPARATraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDPARATraj_Solver::*)(const int&)>(&PIMDPARATraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDPARATraj_Solver::*)(const int&)>(&PIMDPARATraj_Solver::rst_read))
    .def("mpiSendx", [](PIMDPARATraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> x_mpi_arr) {
            return self.mpiSendx(x_mpi_arr.mutable_data()); 
        }
    )
    .def("mpiRecvf", [](PIMDPARATraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> f_mpi_arr, 
        int nsize) {
            return self.mpiRecvf(f_mpi_arr.mutable_data(), nsize); 
        }
    )
    .def("printdata", static_cast<int (PIMDPARATraj_Solver::*)()>(&PIMDPARATraj_Solver::printdata))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_MBPIMDTraj_Solver : public MBPIMDTraj_Solver {
        public:
        using MBPIMDTraj_Solver::MBPIMDTraj_Solver;

        int spring_force() override {
            PYBIND11_OVERRIDE(
            int, // return type
            MBPIMDTraj_Solver, // parent class
            spring_force, // func name
            
            );
        }

        int update_p_harm(const num_real &dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MBPIMDTraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int estimator(const int &isamp, TCFnucl &tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MBPIMDTraj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<MBPIMDTraj_Solver, PIMDTraj_Solver, PyTrampoline_MBPIMDTraj_Solver>(solvers_m, "MBPIMDTraj_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_fV", &MBPIMDTraj_Solver::ref_fV, py::return_value_policy::reference_internal)
    .def("ref_fE", &MBPIMDTraj_Solver::ref_fE, py::return_value_policy::reference_internal)
    .def("ref_VHO", &MBPIMDTraj_Solver::ref_VHO, py::return_value_policy::reference_internal)
    .def("ref_dV_spring", &MBPIMDTraj_Solver::ref_dV_spring, py::return_value_policy::reference_internal)
    .def("ref_dE_spring", &MBPIMDTraj_Solver::ref_dE_spring, py::return_value_policy::reference_internal)
    .def_static("name", &MBPIMDTraj_Solver::name)
    .def("spring_force", static_cast<int (MBPIMDTraj_Solver::*)()>(&MBPIMDTraj_Solver::spring_force))
    .def("update_p_harm", static_cast<int (MBPIMDTraj_Solver::*)(const num_real&)>(&MBPIMDTraj_Solver::update_p_harm))
    .def("estimator", static_cast<int (MBPIMDTraj_Solver::*)(const int&, TCFnucl&)>(&MBPIMDTraj_Solver::estimator))
    .def("ff_calc1", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::init))
    .def("all_X2K", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDTraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("pseudo_inv", [](PIMDTraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_read))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_PILD_Solver : public PILD_Solver {
        public:
        using PILD_Solver::PILD_Solver;

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PILD_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PILD_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PILD_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PILD_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p_harm(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<PILD_Solver, PIMDTraj_Solver, PyTrampoline_PILD_Solver>(solvers_m, "PILD_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_M_therm", &PILD_Solver::ref_M_therm, py::return_value_policy::reference_internal)
    .def_static("name", &PILD_Solver::name)
    .def("update_thermo", static_cast<int (PILD_Solver::*)(const num_real&)>(&PILD_Solver::update_thermo))
    .def("update_p", static_cast<int (PILD_Solver::*)(const num_real&)>(&PILD_Solver::update_p))
    .def("init", static_cast<int (PILD_Solver::*)(const int&)>(&PILD_Solver::init))
    .def("estimator", static_cast<int (PILD_Solver::*)(const int&, TCFnucl&)>(&PILD_Solver::estimator))
    .def("ff_calc1", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::init))
    .def("all_X2K", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDTraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("pseudo_inv", [](PIMDTraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_read))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_CentroidMD_Solver : public CentroidMD_Solver {
        public:
        using CentroidMD_Solver::CentroidMD_Solver;

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CentroidMD_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CentroidMD_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CentroidMD_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_p_harm(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            update_p_harm, // func name
            dt
            );
        }

        int caylay_update_half(
        const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            caylay_update_half, // func name
            dt
            );
        }

        int BAOAB(int& succ,
                      int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BAOAB, // func name
            succ, step
            );
        }

        int BCOCB(int& succ, int step) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            BCOCB, // func name
            succ, step
            );
        }

        int traj(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj, // func name
            tcfer, PN
            );
        }

        int traj_Middle(TCFnucl& tcfer, const int& PN) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_Middle, // func name
            tcfer, PN
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rot_trans_corr, // func name
            Natom, m_in, x_in, p_in, F_in, cal_force
            );
        }

        int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            pseudo_inv, // func name
            N, A, invA, vectmp, eps
            );
        }

        int cross(num_real* vec1, num_real* vec2, num_real* prod) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            cross, // func name
            vec1, vec2, prod
            );
        }

        int rst_output(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_output, // func name
            traj_in
            );
        }

        int rst_read(const int& traj_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PIMDTraj_Solver, // parent class
            rst_read, // func name
            traj_in
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<CentroidMD_Solver, PIMDTraj_Solver, PyTrampoline_CentroidMD_Solver>(solvers_m, "CentroidMD_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &CentroidMD_Solver::name)
    .def("update_thermo", static_cast<int (CentroidMD_Solver::*)(const num_real&)>(&CentroidMD_Solver::update_thermo))
    .def("init", static_cast<int (CentroidMD_Solver::*)(const int&)>(&CentroidMD_Solver::init))
    .def("estimator", static_cast<int (CentroidMD_Solver::*)(const int&, TCFnucl&)>(&CentroidMD_Solver::estimator))
    .def("ff_calc1", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::ff_calc1))
    .def("update_r", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_r))
    .def("update_p", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p))
    .def("update_p_harm", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_p_harm))
    .def("caylay_update_half", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::caylay_update_half))
    .def("BAOAB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BAOAB))
    .def("BCOCB", static_cast<int (PIMDTraj_Solver::*)(int&, int)>(&PIMDTraj_Solver::BCOCB))
    .def("update_thermo", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::update_thermo))
    .def("traj", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj))
    .def("traj_Middle", static_cast<int (PIMDTraj_Solver::*)(TCFnucl&, const int&)>(&PIMDTraj_Solver::traj_Middle))
    .def("traj_property", static_cast<int (PIMDTraj_Solver::*)(const num_real&)>(&PIMDTraj_Solver::traj_property))
    .def("sampler", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::sampler))
    .def("estimator", static_cast<int (PIMDTraj_Solver::*)(const int&, TCFnucl&)>(&PIMDTraj_Solver::estimator))
    .def("run_impl", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::run_parallel))
    .def("init", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::init))
    .def("all_X2K", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_X2K))
    .def("all_K2X", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_K2X))
    .def("all_FX2FK", static_cast<int (PIMDTraj_Solver::*)()>(&PIMDTraj_Solver::all_FX2FK))
    .def("rot_trans_corr", [](PIMDTraj_Solver& self, int Natom, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> m_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> x_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> p_in_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> F_in_arr, 
        bool cal_force) {
            return self.rot_trans_corr(Natom, m_in_arr.mutable_data(), x_in_arr.mutable_data(), p_in_arr.mutable_data(), F_in_arr.mutable_data(), cal_force); 
        }
    )
    .def("pseudo_inv", [](PIMDTraj_Solver& self, int N, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> A_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> invA_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vectmp_arr, 
        num_real eps) {
            return self.pseudo_inv(N, A_arr.mutable_data(), invA_arr.mutable_data(), vectmp_arr.mutable_data(), eps); 
        }
    )
    .def("cross", [](PIMDTraj_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> vec1_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> vec2_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> prod_arr) {
            return self.cross(vec1_arr.mutable_data(), vec2_arr.mutable_data(), prod_arr.mutable_data()); 
        }
    )
    .def("rst_output", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (PIMDTraj_Solver::*)(const int&)>(&PIMDTraj_Solver::rst_read))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_NadTraj_Solver : public NadTraj_Solver {
        public:
        using NadTraj_Solver::NadTraj_Solver;

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<NadTraj_Solver, Traj_Solver, PyTrampoline_NadTraj_Solver>(solvers_m, "NadTraj_Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_mvc", &NadTraj_Solver::ref_mvc, py::return_value_policy::reference_internal)
    .def("ref_fmean", &NadTraj_Solver::ref_fmean, py::return_value_policy::reference_internal)
    .def("ref_fcorr", &NadTraj_Solver::ref_fcorr, py::return_value_policy::reference_internal)
    .def("ref_V", &NadTraj_Solver::ref_V, py::return_value_policy::reference_internal)
    .def("ref_dV", &NadTraj_Solver::ref_dV, py::return_value_policy::reference_internal)
    .def("ref_ddV", &NadTraj_Solver::ref_ddV, py::return_value_policy::reference_internal)
    .def("ref_E", &NadTraj_Solver::ref_E, py::return_value_policy::reference_internal)
    .def("ref_T", &NadTraj_Solver::ref_T, py::return_value_policy::reference_internal)
    .def("ref_T0", &NadTraj_Solver::ref_T0, py::return_value_policy::reference_internal)
    .def("ref_dE", &NadTraj_Solver::ref_dE, py::return_value_policy::reference_internal)
    .def("ref_ddE", &NadTraj_Solver::ref_ddE, py::return_value_policy::reference_internal)
    .def("ref_L", &NadTraj_Solver::ref_L, py::return_value_policy::reference_internal)
    .def("ref_nacv", &NadTraj_Solver::ref_nacv, py::return_value_policy::reference_internal)
    .def("ref_eac0", &NadTraj_Solver::ref_eac0, py::return_value_policy::reference_internal)
    .def("ref_rho0", &NadTraj_Solver::ref_rho0, py::return_value_policy::reference_internal)
    .def("ref_rhot", &NadTraj_Solver::ref_rhot, py::return_value_policy::reference_internal)
    .def("ref_U", &NadTraj_Solver::ref_U, py::return_value_policy::reference_internal)
    .def("ref_H", &NadTraj_Solver::ref_H, py::return_value_policy::reference_internal)
    .def("ref_dH", &NadTraj_Solver::ref_dH, py::return_value_policy::reference_internal)
    .def("ref_ddH", &NadTraj_Solver::ref_ddH, py::return_value_policy::reference_internal)
    .def("ref_S", &NadTraj_Solver::ref_S, py::return_value_policy::reference_internal)
    .def("ref_dL", &NadTraj_Solver::ref_dL, py::return_value_policy::reference_internal)
    .def("ref_ddL", &NadTraj_Solver::ref_ddL, py::return_value_policy::reference_internal)
    .def("ref_eac", &NadTraj_Solver::ref_eac, py::return_value_policy::reference_internal)
    .def("ref_eacf", &NadTraj_Solver::ref_eacf, py::return_value_policy::reference_internal)
    .def("ref_eacb", &NadTraj_Solver::ref_eacb, py::return_value_policy::reference_internal)
    .def("ref_gmat", &NadTraj_Solver::ref_gmat, py::return_value_policy::reference_internal)
    .def("ref_rho", &NadTraj_Solver::ref_rho, py::return_value_policy::reference_internal)
    .def_static("name", &NadTraj_Solver::name)
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_MMD_Solver : public MMD_Solver {
        public:
        using MMD_Solver::MMD_Solver;

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MMD_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MMD_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            MMD_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<MMD_Solver, NadTraj_Solver, PyTrampoline_MMD_Solver>(solvers_m, "MMD_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &MMD_Solver::name)
    .def("init", static_cast<int (MMD_Solver::*)(const int&)>(&MMD_Solver::init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_CMM_Solver : public CMM_Solver {
        public:
        using CMM_Solver::CMM_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CMM_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CMM_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CMM_Solver, // parent class
            kernel_cmm, // func name
            rhox, xic, gammac, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CMM_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            CMM_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<CMM_Solver, NadTraj_Solver, PyTrampoline_CMM_Solver>(solvers_m, "CMM_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &CMM_Solver::name)
    .def("init_occ2eac", static_cast<int (CMM_Solver::*)(const int&)>(&CMM_Solver::init_occ2eac))
    .def("init", static_cast<int (CMM_Solver::*)(const int&)>(&CMM_Solver::init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_SQC_Solver : public SQC_Solver {
        public:
        using SQC_Solver::SQC_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SQC_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SQC_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SQC_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SQC_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<SQC_Solver, NadTraj_Solver, PyTrampoline_SQC_Solver>(solvers_m, "SQC_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &SQC_Solver::name)
    .def("init_occ2eac", static_cast<int (SQC_Solver::*)(const int&)>(&SQC_Solver::init_occ2eac))
    .def("init", static_cast<int (SQC_Solver::*)(const int&)>(&SQC_Solver::init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_LSC_Solver : public LSC_Solver {
        public:
        using LSC_Solver::LSC_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LSC_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LSC_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel_lsc(num_complex* rhox, num_real& gammac, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LSC_Solver, // parent class
            kernel_lsc, // func name
            rhox, gammac, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LSC_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            LSC_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<LSC_Solver, NadTraj_Solver, PyTrampoline_LSC_Solver>(solvers_m, "LSC_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &LSC_Solver::name)
    .def("init_occ2eac", static_cast<int (LSC_Solver::*)(const int&)>(&LSC_Solver::init_occ2eac))
    .def("init", static_cast<int (LSC_Solver::*)(const int&)>(&LSC_Solver::init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_SMM_Solver : public SMM_Solver {
        public:
        using SMM_Solver::SMM_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SMM_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SMM_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel_scmm(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SMM_Solver, // parent class
            kernel_scmm, // func name
            rhox, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SMM_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SMM_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<SMM_Solver, NadTraj_Solver, PyTrampoline_SMM_Solver>(solvers_m, "SMM_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_xi1", &SMM_Solver::ref_xi1, py::return_value_policy::reference_internal)
    .def_static("name", &SMM_Solver::name)
    .def("init_occ2eac", static_cast<int (SMM_Solver::*)(const int&)>(&SMM_Solver::init_occ2eac))
    .def("init", static_cast<int (SMM_Solver::*)(const int&)>(&SMM_Solver::init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_Multi_NadTraj_Solver : public Multi_NadTraj_Solver {
        public:
        using Multi_NadTraj_Solver::Multi_NadTraj_Solver;

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int multi_ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            multi_ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int multi_ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            multi_ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int multi_evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            multi_evolve_elec, // func name
            Uevolve
            );
        }

        int update_p(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_p, // func name
            dt_in
            );
        }

        int update_r(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_r, // func name
            dt_in
            );
        }

        int update_thermo(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_thermo, // func name
            dt_in
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int multi_init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            multi_init, // func name
            itraj
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<Multi_NadTraj_Solver, NadTraj_Solver, PyTrampoline_Multi_NadTraj_Solver>(solvers_m, "Multi_NadTraj_Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*, const int&>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_nrs", &Multi_NadTraj_Solver::ref_nrs, py::return_value_policy::reference_internal)
    .def("ref_nps", &Multi_NadTraj_Solver::ref_nps, py::return_value_policy::reference_internal)
    .def("ref_nms", &Multi_NadTraj_Solver::ref_nms, py::return_value_policy::reference_internal)
    .def("ref_nfs", &Multi_NadTraj_Solver::ref_nfs, py::return_value_policy::reference_internal)
    .def("ref_vpeses", &Multi_NadTraj_Solver::ref_vpeses, py::return_value_policy::reference_internal)
    .def("ref_grads", &Multi_NadTraj_Solver::ref_grads, py::return_value_policy::reference_internal)
    .def("ref_hesses", &Multi_NadTraj_Solver::ref_hesses, py::return_value_policy::reference_internal)
    .def("ref_Vs", &Multi_NadTraj_Solver::ref_Vs, py::return_value_policy::reference_internal)
    .def("ref_dVs", &Multi_NadTraj_Solver::ref_dVs, py::return_value_policy::reference_internal)
    .def("ref_ddVs", &Multi_NadTraj_Solver::ref_ddVs, py::return_value_policy::reference_internal)
    .def("ref_Es", &Multi_NadTraj_Solver::ref_Es, py::return_value_policy::reference_internal)
    .def("ref_Ts", &Multi_NadTraj_Solver::ref_Ts, py::return_value_policy::reference_internal)
    .def("ref_dEs", &Multi_NadTraj_Solver::ref_dEs, py::return_value_policy::reference_internal)
    .def("ref_nacvs", &Multi_NadTraj_Solver::ref_nacvs, py::return_value_policy::reference_internal)
    .def("ref_Ls", &Multi_NadTraj_Solver::ref_Ls, py::return_value_policy::reference_internal)
    .def("ref_Ss", &Multi_NadTraj_Solver::ref_Ss, py::return_value_policy::reference_internal)
    .def("ref_rhos", &Multi_NadTraj_Solver::ref_rhos, py::return_value_policy::reference_internal)
    .def("ref_Us", &Multi_NadTraj_Solver::ref_Us, py::return_value_policy::reference_internal)
    .def_static("name", &Multi_NadTraj_Solver::name)
    .def("ff_calc1", static_cast<int (Multi_NadTraj_Solver::*)(const int&, const bool&)>(&Multi_NadTraj_Solver::ff_calc1))
    .def("multi_ff_calc1", static_cast<int (Multi_NadTraj_Solver::*)(const int&, const bool&)>(&Multi_NadTraj_Solver::multi_ff_calc1))
    .def("ff_calc2", static_cast<int (Multi_NadTraj_Solver::*)()>(&Multi_NadTraj_Solver::ff_calc2))
    .def("multi_ff_calc2", static_cast<int (Multi_NadTraj_Solver::*)()>(&Multi_NadTraj_Solver::multi_ff_calc2))
    .def("update_p", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_p))
    .def("update_r", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_r))
    .def("update_thermo", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_thermo))
    .def("init", static_cast<int (Multi_NadTraj_Solver::*)(const int&)>(&Multi_NadTraj_Solver::init))
    .def("multi_init", static_cast<int (Multi_NadTraj_Solver::*)(const int&)>(&Multi_NadTraj_Solver::multi_init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_PMM_Solver : public PMM_Solver {
        public:
        using PMM_Solver::PMM_Solver;

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int multi_evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            multi_evolve_elec, // func name
            Uevolve
            );
        }

        int multi_ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            multi_ff_calc2, // func name
            
            );
        }

        int multi_init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            multi_init, // func name
            itraj
            );
        }

        int multi_reinit(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            multi_reinit, // func name
            itraj
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            PMM_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int multi_ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            multi_ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int update_p(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_p, // func name
            dt_in
            );
        }

        int update_r(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_r, // func name
            dt_in
            );
        }

        int update_thermo(const num_real& dt_in) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            update_thermo, // func name
            dt_in
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Multi_NadTraj_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<PMM_Solver, Multi_NadTraj_Solver, PyTrampoline_PMM_Solver>(solvers_m, "PMM_Solver", py::dynamic_attr())
    .def(py::init<const Param&, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_drhos", &PMM_Solver::ref_drhos, py::return_value_policy::reference_internal)
    .def("ref_deltarhos", &PMM_Solver::ref_deltarhos, py::return_value_policy::reference_internal)
    .def("ref_rhosumt", &PMM_Solver::ref_rhosumt, py::return_value_policy::reference_internal)
    .def("ref_rhoavg", &PMM_Solver::ref_rhoavg, py::return_value_policy::reference_internal)
    .def("ref_rho_corr", &PMM_Solver::ref_rho_corr, py::return_value_policy::reference_internal)
    .def_static("name", &PMM_Solver::name)
    .def("traj", static_cast<int (PMM_Solver::*)(NAD_TCFer&, const int&, const int&)>(&PMM_Solver::traj))
    .def("multi_ff_calc2", static_cast<int (PMM_Solver::*)()>(&PMM_Solver::multi_ff_calc2))
    .def("multi_init", static_cast<int (PMM_Solver::*)(const int&)>(&PMM_Solver::multi_init))
    .def("multi_reinit", static_cast<int (PMM_Solver::*)(const int&)>(&PMM_Solver::multi_reinit))
    .def("ff_calc1", static_cast<int (Multi_NadTraj_Solver::*)(const int&, const bool&)>(&Multi_NadTraj_Solver::ff_calc1))
    .def("multi_ff_calc1", static_cast<int (Multi_NadTraj_Solver::*)(const int&, const bool&)>(&Multi_NadTraj_Solver::multi_ff_calc1))
    .def("ff_calc2", static_cast<int (Multi_NadTraj_Solver::*)()>(&Multi_NadTraj_Solver::ff_calc2))
    .def("multi_ff_calc2", static_cast<int (Multi_NadTraj_Solver::*)()>(&Multi_NadTraj_Solver::multi_ff_calc2))
    .def("update_p", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_p))
    .def("update_r", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_r))
    .def("update_thermo", static_cast<int (Multi_NadTraj_Solver::*)(const num_real&)>(&Multi_NadTraj_Solver::update_thermo))
    .def("init", static_cast<int (Multi_NadTraj_Solver::*)(const int&)>(&Multi_NadTraj_Solver::init))
    .def("multi_init", static_cast<int (Multi_NadTraj_Solver::*)(const int&)>(&Multi_NadTraj_Solver::multi_init))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_WMM_Solver : public WMM_Solver {
        public:
        using WMM_Solver::WMM_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel_cmm(num_complex* rhox, num_real& xic, num_real& gammac, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            kernel_cmm, // func name
            rhox, xic, gammac, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            WMM_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<WMM_Solver, NadTraj_Solver, PyTrampoline_WMM_Solver>(solvers_m, "WMM_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_static("name", &WMM_Solver::name)
    .def("init_occ2eac", static_cast<int (WMM_Solver::*)(const int&)>(&WMM_Solver::init_occ2eac))
    .def("init", static_cast<int (WMM_Solver::*)(const int&)>(&WMM_Solver::init))
    .def("traj_velocityverlet", static_cast<int (WMM_Solver::*)(NAD_TCFer&, const int&, const int&)>(&WMM_Solver::traj_velocityverlet))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_NAD_TCFer : public NAD_TCFer {
        public:
        using NAD_TCFer::NAD_TCFer;
    };

    py::class_<NAD_TCFer, PyTrampoline_NAD_TCFer>(models_m, "NAD_TCFer", py::dynamic_attr())
    .def(py::init<const int&, const int&, const int&, const int&, const int&, const int&>())
    .def_readwrite("lentcf", &NAD_TCFer::lentcf)
    .def_readwrite("tcf_reduced", &NAD_TCFer::tcf_reduced)
    .def("ref_val", &NAD_TCFer::ref_val, py::return_value_policy::reference_internal)
    .def("ref_tcf_0t_bool", &NAD_TCFer::ref_tcf_0t_bool, py::return_value_policy::reference_internal)
    .def("ref_tcf_0t_val", &NAD_TCFer::ref_tcf_0t_val, py::return_value_policy::reference_internal)
    .def("ref_coll", &NAD_TCFer::ref_coll, py::return_value_policy::reference_internal)
    .def("ref_stat", &NAD_TCFer::ref_stat, py::return_value_policy::reference_internal)
    .def("ref_tcf_0_bool", &NAD_TCFer::ref_tcf_0_bool, py::return_value_policy::reference_internal)
    .def("ref_tcf_t_bool", &NAD_TCFer::ref_tcf_t_bool, py::return_value_policy::reference_internal)
    .def("ref_tcf_0_val", &NAD_TCFer::ref_tcf_0_val, py::return_value_policy::reference_internal)
    .def("ref_tcf_t_val", &NAD_TCFer::ref_tcf_t_val, py::return_value_policy::reference_internal)
    .def("Clear", static_cast<int (NAD_TCFer::*)()>(&NAD_TCFer::Clear))
    .def("Amount", static_cast<int (NAD_TCFer::*)(NAD_TCFer&)>(&NAD_TCFer::Amount))
    .def("MPIAmount", static_cast<int (NAD_TCFer::*)(NAD_TCFer&)>(&NAD_TCFer::MPIAmount))
    .def("Count", static_cast<int (NAD_TCFer::*)(const int&, const int&)>(&NAD_TCFer::Count))
    .def("ifrecord", static_cast<bool (NAD_TCFer::*)(const int&, const int&)>(&NAD_TCFer::ifrecord))
    .def("report", static_cast<int (NAD_TCFer::*)(const std::string&, const double&)>(&NAD_TCFer::report));

    class PyTrampoline_Hopping_Solver : public Hopping_Solver {
        public:
        using Hopping_Solver::Hopping_Solver;

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            init, // func name
            itraj
            );
        }

        inline num_real nucl_Ekin() override {
            PYBIND11_OVERRIDE(
            num_real, // return type
            Hopping_Solver, // parent class
            nucl_Ekin, // func name
            
            );
        }

        int SE_Hamiltonian(const bool& refered) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            SE_Hamiltonian, // func name
            refered
            );
        }

        int phase_correction() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            phase_correction, // func name
            
            );
        }

        int hopping_choose(const int& iocc, num_complex* rhox, num_complex* H, num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            hopping_choose, // func name
            iocc, rhox, H, dt
            );
        }

        int hopping_direction(num_real* direction, const int& to) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            hopping_direction, // func name
            direction, to
            );
        }

        int hopping_impulse(num_real* direction, const int& occ_to) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            hopping_impulse, // func name
            direction, occ_to
            );
        }

        int time_calc() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            time_calc, // func name
            
            );
        }

        int coherent_evolve(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            coherent_evolve, // func name
            dt
            );
        }

        int decoherent_evolve(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            decoherent_evolve, // func name
            dt
            );
        }

        int hopping_collapse(const num_real& dt, const int& to) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            hopping_collapse, // func name
            dt, to
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int kernel_fssh(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            kernel_fssh, // func name
            rhox, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Hopping_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<Hopping_Solver, NadTraj_Solver, PyTrampoline_Hopping_Solver>(solvers_m, "Hopping_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_arg_nr", &Hopping_Solver::ref_arg_nr, py::return_value_policy::reference_internal)
    .def("ref_arg_np", &Hopping_Solver::ref_arg_np, py::return_value_policy::reference_internal)
    .def("ref_time_coh", &Hopping_Solver::ref_time_coh, py::return_value_policy::reference_internal)
    .def("ref_time_los", &Hopping_Solver::ref_time_los, py::return_value_policy::reference_internal)
    .def("ref_dnp_diss", &Hopping_Solver::ref_dnp_diss, py::return_value_policy::reference_internal)
    .def("ref_probs_dish", &Hopping_Solver::ref_probs_dish, py::return_value_policy::reference_internal)
    .def("ref_idx_dish", &Hopping_Solver::ref_idx_dish, py::return_value_policy::reference_internal)
    .def("ref_direction", &Hopping_Solver::ref_direction, py::return_value_policy::reference_internal)
    .def("ref_drho_diss", &Hopping_Solver::ref_drho_diss, py::return_value_policy::reference_internal)
    .def("ref_Uh", &Hopping_Solver::ref_Uh, py::return_value_policy::reference_internal)
    .def_static("name", &Hopping_Solver::name)
    .def("init_occ2eac", static_cast<int (Hopping_Solver::*)(const int&)>(&Hopping_Solver::init_occ2eac))
    .def("init", static_cast<int (Hopping_Solver::*)(const int&)>(&Hopping_Solver::init))
    .def("nucl_Ekin", static_cast<num_real (Hopping_Solver::*)()>(&Hopping_Solver::nucl_Ekin))
    .def("SE_Hamiltonian", static_cast<int (Hopping_Solver::*)(const bool&)>(&Hopping_Solver::SE_Hamiltonian))
    .def("phase_correction", static_cast<int (Hopping_Solver::*)()>(&Hopping_Solver::phase_correction))
    .def("hopping_direction", [](Hopping_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> direction_arr, 
        const int& to) {
            return self.hopping_direction(direction_arr.mutable_data(), to); 
        }
    )
    .def("hopping_impulse", [](Hopping_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> direction_arr, 
        const int& occ_to) {
            return self.hopping_impulse(direction_arr.mutable_data(), occ_to); 
        }
    )
    .def("time_calc", static_cast<int (Hopping_Solver::*)()>(&Hopping_Solver::time_calc))
    .def("coherent_evolve", static_cast<int (Hopping_Solver::*)(const num_real&)>(&Hopping_Solver::coherent_evolve))
    .def("decoherent_evolve", static_cast<int (Hopping_Solver::*)(const num_real&)>(&Hopping_Solver::decoherent_evolve))
    .def("hopping_collapse", static_cast<int (Hopping_Solver::*)(const num_real&, const int&)>(&Hopping_Solver::hopping_collapse))
    .def("ff_calc1", static_cast<int (Hopping_Solver::*)(const int&, const bool&)>(&Hopping_Solver::ff_calc1))
    .def("check_break", static_cast<int (Hopping_Solver::*)(int&)>(&Hopping_Solver::check_break))
    .def("ff_calc2", static_cast<int (Hopping_Solver::*)()>(&Hopping_Solver::ff_calc2))
    .def("traj", static_cast<int (Hopping_Solver::*)(NAD_TCFer&, const int&, const int&)>(&Hopping_Solver::traj))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_ProductMPS_Solver : public ProductMPS_Solver {
        public:
        using ProductMPS_Solver::ProductMPS_Solver;

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ProductMPS_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ProductMPS_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            ProductMPS_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<ProductMPS_Solver, NadTraj_Solver, PyTrampoline_ProductMPS_Solver>(solvers_m, "ProductMPS_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_readwrite("M", &ProductMPS_Solver::M)
    .def_readwrite("F", &ProductMPS_Solver::F)
    .def_readwrite("FF", &ProductMPS_Solver::FF)
    .def_readwrite("MF", &ProductMPS_Solver::MF)
    .def_readwrite("MFF", &ProductMPS_Solver::MFF)
    .def_readwrite("type", &ProductMPS_Solver::type)
    .def_readwrite("scale", &ProductMPS_Solver::scale)
    .def("ref_rhos", &ProductMPS_Solver::ref_rhos, py::return_value_policy::reference_internal)
    .def("ref_rhos0", &ProductMPS_Solver::ref_rhos0, py::return_value_policy::reference_internal)
    .def("ref_rhost", &ProductMPS_Solver::ref_rhost, py::return_value_policy::reference_internal)
    .def("ref_Hs", &ProductMPS_Solver::ref_Hs, py::return_value_policy::reference_internal)
    .def_static("name", &ProductMPS_Solver::name)
    .def("init", static_cast<int (ProductMPS_Solver::*)(const int&)>(&ProductMPS_Solver::init))
    .def("traj", static_cast<int (ProductMPS_Solver::*)(NAD_TCFer&, const int&, const int&)>(&ProductMPS_Solver::traj))
    .def("correlation", static_cast<int (ProductMPS_Solver::*)(const int&, NAD_TCFer&)>(&ProductMPS_Solver::correlation))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_QCPI_Solver : public QCPI_Solver {
        public:
        using QCPI_Solver::QCPI_Solver;

        int propagate_tensor() override {
            PYBIND11_OVERRIDE(
            int, // return type
            QCPI_Solver, // parent class
            propagate_tensor, // func name
            
            );
        }

        num_complex propagate_nuc(num_real* nr2_trj, num_real* np2_trj,  
                                      num_real* nr1_trj, num_real* np1_trj,  
                                      const int& prev, const int& next, const num_real& dt) override {
            PYBIND11_OVERRIDE(
            num_complex, // return type
            QCPI_Solver, // parent class
            propagate_nuc, // func name
            nr2_trj, np2_trj, nr1_trj, np1_trj, prev, next, dt
            );
        }

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            QCPI_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int kernel0(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            QCPI_Solver, // parent class
            kernel0, // func name
            rhox, F
            );
        }

        int kernelt(num_complex* rhox, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            QCPI_Solver, // parent class
            kernelt, // func name
            rhox, F
            );
        }

        int traj(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            QCPI_Solver, // parent class
            traj, // func name
            tcfer, N, F
            );
        }

        int ff_calc1(const int& level, const bool& refered = false) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc1, // func name
            level, refered
            );
        }

        int ff_calc2() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            ff_calc2, // func name
            
            );
        }

        int evolve_elec(num_complex* Uevolve) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            evolve_elec, // func name
            Uevolve
            );
        }

        int init_occ2eac(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_occ2eac, // func name
            itraj
            );
        }

        int init_ofs(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            init_ofs, // func name
            itraj
            );
        }

        int final(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            final, // func name
            itraj
            );
        }

        int check_break(int& succ) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            check_break, // func name
            succ
            );
        }

        int rst_output(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_output, // func name
            itraj
            );
        }

        int rst_read(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            rst_read, // func name
            itraj
            );
        }

        int traj_property(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_property, // func name
            dt
            );
        }

        int traj_velocityverlet(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N, F
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            NadTraj_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int ff_calc1(const int& level) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            ff_calc1, // func name
            level
            );
        }

        int update_r(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_r, // func name
            dt
            );
        }

        int update_p(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_p, // func name
            dt
            );
        }

        int update_thermo(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            update_thermo, // func name
            dt
            );
        }

        int traj(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj, // func name
            tcfer, N
            );
        }

        int traj_velocityverlet(TCFnucl& tcfer, const int& N) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            traj_velocityverlet, // func name
            tcfer, N
            );
        }

        int sampler(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int estimator(const int& isamp, TCFnucl& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Traj_Solver, // parent class
            estimator, // func name
            isamp, tcfer
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<QCPI_Solver, NadTraj_Solver, PyTrampoline_QCPI_Solver>(solvers_m, "QCPI_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_idx_arr", &QCPI_Solver::ref_idx_arr, py::return_value_policy::reference_internal)
    .def("ref_nrs", &QCPI_Solver::ref_nrs, py::return_value_policy::reference_internal)
    .def("ref_nps", &QCPI_Solver::ref_nps, py::return_value_policy::reference_internal)
    .def("ref_nrs_copy", &QCPI_Solver::ref_nrs_copy, py::return_value_policy::reference_internal)
    .def("ref_nps_copy", &QCPI_Solver::ref_nps_copy, py::return_value_policy::reference_internal)
    .def("ref_nx", &QCPI_Solver::ref_nx, py::return_value_policy::reference_internal)
    .def("ref_ny", &QCPI_Solver::ref_ny, py::return_value_policy::reference_internal)
    .def("ref_nx_copy", &QCPI_Solver::ref_nx_copy, py::return_value_policy::reference_internal)
    .def("ref_ny_copy", &QCPI_Solver::ref_ny_copy, py::return_value_policy::reference_internal)
    .def_static("name", &QCPI_Solver::name)
    .def("get_index", [](QCPI_Solver& self, py::array_t<int, py::array::c_style | py::array::forcecast> arr_arr, 
        const int& num, 
        const int& base, 
        const int& len) {
            return self.get_index(arr_arr.mutable_data(), num, base, len); 
        }
    )
    .def("get_Nnum", [](QCPI_Solver& self, py::array_t<int, py::array::c_style | py::array::forcecast> arr_arr, 
        const int& base, 
        const int& len) {
            return self.get_Nnum(arr_arr.mutable_data(), base, len); 
        }
    )
    .def("propagate_tensor", static_cast<int (QCPI_Solver::*)()>(&QCPI_Solver::propagate_tensor))
    .def("propagate_nuc", [](QCPI_Solver& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr2_trj_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np2_trj_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nr1_trj_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np1_trj_arr, 
        const int& prev, 
        const int& next, 
        const num_real& dt) {
            return self.propagate_nuc(nr2_trj_arr.mutable_data(), np2_trj_arr.mutable_data(), nr1_trj_arr.mutable_data(), np1_trj_arr.mutable_data(), prev, next, dt); 
        }
    )
    .def("init", static_cast<int (QCPI_Solver::*)(const int&)>(&QCPI_Solver::init))
    .def("traj", static_cast<int (QCPI_Solver::*)(NAD_TCFer&, const int&, const int&)>(&QCPI_Solver::traj))
    .def("ff_calc1", static_cast<int (NadTraj_Solver::*)(const int&, const bool&)>(&NadTraj_Solver::ff_calc1))
    .def("ff_calc2", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::ff_calc2))
    .def("init_occ2eac", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_occ2eac))
    .def("init_ofs", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init_ofs))
    .def("init", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::init))
    .def("final", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::final))
    .def("check_break", static_cast<int (NadTraj_Solver::*)(int&)>(&NadTraj_Solver::check_break))
    .def("rst_output", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_output))
    .def("rst_read", static_cast<int (NadTraj_Solver::*)(const int&)>(&NadTraj_Solver::rst_read))
    .def("traj_property", static_cast<int (NadTraj_Solver::*)(const num_real&)>(&NadTraj_Solver::traj_property))
    .def("traj", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (NadTraj_Solver::*)(NAD_TCFer&, const int&, const int&)>(&NadTraj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::sampler))
    .def("correlation", static_cast<int (NadTraj_Solver::*)(const int&, NAD_TCFer&)>(&NadTraj_Solver::correlation))
    .def("run_impl", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_impl))
    .def("run_parallel", static_cast<int (NadTraj_Solver::*)()>(&NadTraj_Solver::run_parallel))
    .def("ff_calc1", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::ff_calc1))
    .def("init_ofs", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init_ofs))
    .def("init", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::init))
    .def("final", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::final))
    .def("rst_output", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_output))
    .def("rst_read", static_cast<int (Traj_Solver::*)(const int&)>(&Traj_Solver::rst_read))
    .def("check_break", static_cast<int (Traj_Solver::*)(int&)>(&Traj_Solver::check_break))
    .def("traj_property", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::traj_property))
    .def("update_r", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_r))
    .def("update_p", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_p))
    .def("update_thermo", static_cast<int (Traj_Solver::*)(const num_real&)>(&Traj_Solver::update_thermo))
    .def("traj", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj))
    .def("traj_velocityverlet", static_cast<int (Traj_Solver::*)(TCFnucl&, const int&)>(&Traj_Solver::traj_velocityverlet))
    .def("sampler", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::sampler))
    .def("estimator", static_cast<int (Traj_Solver::*)(const int&, TCFnucl&)>(&Traj_Solver::estimator))
    .def("run_impl", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_impl))
    .def("run_parallel", static_cast<int (Traj_Solver::*)()>(&Traj_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_HEOM_Solver : public HEOM_Solver {
        public:
        using HEOM_Solver::HEOM_Solver;

        int Init_Bath(const int& ib, const int& Nr, const int& Ni, num_real& gdph_m, num_complex* tcf_coef,
                          num_complex* tcf_zero, num_complex* tcf_deri) override {
            PYBIND11_OVERRIDE(
            int, // return type
            HEOM_Solver, // parent class
            Init_Bath, // func name
            ib, Nr, Ni, gdph_m, tcf_coef, tcf_zero, tcf_deri
            );
        }

        int Generate_ADOs() override {
            PYBIND11_OVERRIDE(
            int, // return type
            HEOM_Solver, // parent class
            Generate_ADOs, // func name
            
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            HEOM_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            HEOM_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int Evolve(num_complex* sigma_tot_t1, num_complex* sigma_tot_t2) override {
            PYBIND11_OVERRIDE(
            int, // return type
            HEOM_Solver, // parent class
            Evolve, // func name
            sigma_tot_t1, sigma_tot_t2
            );
        }

        bool if_filteration(int* arr) override {
            PYBIND11_OVERRIDE(
            bool, // return type
            HEOM_Solver, // parent class
            if_filteration, // func name
            arr
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<HEOM_Solver, Solver, PyTrampoline_HEOM_Solver>(solvers_m, "HEOM_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def_readwrite("N", &HEOM_Solver::N)
    .def_readwrite("H", &HEOM_Solver::H)
    .def_readwrite("F", &HEOM_Solver::F)
    .def_readwrite("FF", &HEOM_Solver::FF)
    .def_readwrite("nbath", &HEOM_Solver::nbath)
    .def_readwrite("Nexpan_M", &HEOM_Solver::Nexpan_M)
    .def_readwrite("Nexpan_Nr", &HEOM_Solver::Nexpan_Nr)
    .def_readwrite("Nexpan_Ni", &HEOM_Solver::Nexpan_Ni)
    .def_readwrite("Nchs", &HEOM_Solver::Nchs)
    .def_readwrite("NADOs", &HEOM_Solver::NADOs)
    .def_readwrite("Nconn", &HEOM_Solver::Nconn)
    .def_readwrite("BasisExpanType", &HEOM_Solver::BasisExpanType)
    .def_readwrite("basis_file", &HEOM_Solver::basis_file)
    .def_readwrite("sstep", &HEOM_Solver::sstep)
    .def_readwrite("nstep", &HEOM_Solver::nstep)
    .def_readwrite("adia_pure", &HEOM_Solver::adia_pure)
    .def_readwrite("occ0", &HEOM_Solver::occ0)
    .def_readwrite("dt", &HEOM_Solver::dt)
    .def_readwrite("tend", &HEOM_Solver::tend)
    .def_readwrite("vscale", &HEOM_Solver::vscale)
    .def("ref_csr_iADOs", &HEOM_Solver::ref_csr_iADOs, py::return_value_policy::reference_internal)
    .def("ref_csr_connec_LD", &HEOM_Solver::ref_csr_connec_LD, py::return_value_policy::reference_internal)
    .def("ref_csr_connec", &HEOM_Solver::ref_csr_connec, py::return_value_policy::reference_internal)
    .def("ref_csr_type", &HEOM_Solver::ref_csr_type, py::return_value_policy::reference_internal)
    .def("ref_csr_ivalue", &HEOM_Solver::ref_csr_ivalue, py::return_value_policy::reference_internal)
    .def("ref_csr_cvalue", &HEOM_Solver::ref_csr_cvalue, py::return_value_policy::reference_internal)
    .def("ref_Tcftype", &HEOM_Solver::ref_Tcftype, py::return_value_policy::reference_internal)
    .def("ref_tcf_site", &HEOM_Solver::ref_tcf_site, py::return_value_policy::reference_internal)
    .def("ref_Esys", &HEOM_Solver::ref_Esys, py::return_value_policy::reference_internal)
    .def("ref_tcf_coef", &HEOM_Solver::ref_tcf_coef, py::return_value_policy::reference_internal)
    .def("ref_tcf_zero", &HEOM_Solver::ref_tcf_zero, py::return_value_policy::reference_internal)
    .def("ref_tcf_deri", &HEOM_Solver::ref_tcf_deri, py::return_value_policy::reference_internal)
    .def("ref_sigma_tot", &HEOM_Solver::ref_sigma_tot, py::return_value_policy::reference_internal)
    .def("ref_rdm_arr", &HEOM_Solver::ref_rdm_arr, py::return_value_policy::reference_internal)
    .def("ref_Tsys", &HEOM_Solver::ref_Tsys, py::return_value_policy::reference_internal)
    .def("ref_Hsys", &HEOM_Solver::ref_Hsys, py::return_value_policy::reference_internal)
    .def("ref_T0", &HEOM_Solver::ref_T0, py::return_value_policy::reference_internal)
    .def("ref_Hsys_time", &HEOM_Solver::ref_Hsys_time, py::return_value_policy::reference_internal)
    .def("ref_eac0", &HEOM_Solver::ref_eac0, py::return_value_policy::reference_internal)
    .def("ref_rho0", &HEOM_Solver::ref_rho0, py::return_value_policy::reference_internal)
    .def_static("name", &HEOM_Solver::name)
    .def("Generate_ADOs", static_cast<int (HEOM_Solver::*)()>(&HEOM_Solver::Generate_ADOs))
    .def("run_impl", static_cast<int (HEOM_Solver::*)()>(&HEOM_Solver::run_impl))
    .def("run_parallel", static_cast<int (HEOM_Solver::*)()>(&HEOM_Solver::run_parallel))
    .def("if_filteration", [](HEOM_Solver& self, py::array_t<int, py::array::c_style | py::array::forcecast> arr_arr) {
            return self.if_filteration(arr_arr.mutable_data()); 
        }
    )
    .def("findarr", [](HEOM_Solver& self, py::array_t<int, py::array::c_style | py::array::forcecast> arr_arr, 
        const long long int& iX, 
        const int& N) {
            return self.findarr(arr_arr.mutable_data(), iX, N); 
        }
    )
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_RedField_Solver : public RedField_Solver {
        public:
        using RedField_Solver::RedField_Solver;

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            RedField_Solver, // parent class
            run_impl, // func name
            
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<RedField_Solver, Solver, PyTrampoline_RedField_Solver>(solvers_m, "RedField_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_Eele", &RedField_Solver::ref_Eele, py::return_value_policy::reference_internal)
    .def("ref_Tele", &RedField_Solver::ref_Tele, py::return_value_policy::reference_internal)
    .def("ref_Hele", &RedField_Solver::ref_Hele, py::return_value_policy::reference_internal)
    .def("ref_mDE", &RedField_Solver::ref_mDE, py::return_value_policy::reference_internal)
    .def("ref_Qtran", &RedField_Solver::ref_Qtran, py::return_value_policy::reference_internal)
    .def("ref_C_mDE", &RedField_Solver::ref_C_mDE, py::return_value_policy::reference_internal)
    .def("ref_GM_tensor", &RedField_Solver::ref_GM_tensor, py::return_value_policy::reference_internal)
    .def("ref_R_tensor", &RedField_Solver::ref_R_tensor, py::return_value_policy::reference_internal)
    .def("ref_eac0", &RedField_Solver::ref_eac0, py::return_value_policy::reference_internal)
    .def("ref_rho0", &RedField_Solver::ref_rho0, py::return_value_policy::reference_internal)
    .def("ref_rhoadia", &RedField_Solver::ref_rhoadia, py::return_value_policy::reference_internal)
    .def("ref_rhodia", &RedField_Solver::ref_rhodia, py::return_value_policy::reference_internal)
    .def("ref_rho1", &RedField_Solver::ref_rho1, py::return_value_policy::reference_internal)
    .def("ref_rho2", &RedField_Solver::ref_rho2, py::return_value_policy::reference_internal)
    .def("ref_rho3", &RedField_Solver::ref_rho3, py::return_value_policy::reference_internal)
    .def("ref_rho4", &RedField_Solver::ref_rho4, py::return_value_policy::reference_internal)
    .def("ref_rhotmp", &RedField_Solver::ref_rhotmp, py::return_value_policy::reference_internal)
    .def_static("name", &RedField_Solver::name)
    .def("run_impl", static_cast<int (RedField_Solver::*)()>(&RedField_Solver::run_impl))
    .def("calc_R_tensor", static_cast<int (RedField_Solver::*)()>(&RedField_Solver::calc_R_tensor))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_SSE_Solver : public SSE_Solver {
        public:
        using SSE_Solver::SSE_Solver;

        int init(const int& itraj) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            init, // func name
            itraj
            );
        }

        int correlation(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            correlation, // func name
            isamp, tcfer
            );
        }

        int sampler(const int& isamp, NAD_TCFer& tcfer) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            sampler, // func name
            isamp, tcfer
            );
        }

        int action_on_wavafunction(num_complex* eacnew, num_complex* eac, const num_real& t) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            action_on_wavafunction, // func name
            eacnew, eac, t
            );
        }

        int action_on_wavafunction_adia(num_complex* eacnew, num_complex* eac, const num_real& t) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            action_on_wavafunction_adia, // func name
            eacnew, eac, t
            );
        }

        int time_kernel(num_complex* Ker, const num_real& t) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            time_kernel, // func name
            Ker, t
            );
        }

        int call_memory_array_adia(const num_real& dt) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            call_memory_array_adia, // func name
            dt
            );
        }

        int sse(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            sse, // func name
            tcfer, N, F
            );
        }

        int sse_adia(NAD_TCFer& tcfer, const int& N, const int& F) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            sse_adia, // func name
            tcfer, N, F
            );
        }

        int action_on_wavafunction_adia_mem(num_complex* eacnew, num_complex* eacold, const int& hstep,
                                                const bool& hupdate) override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            action_on_wavafunction_adia_mem, // func name
            eacnew, eacold, hstep, hupdate
            );
        }

        int run_parallel() override {
            PYBIND11_OVERRIDE(
            int, // return type
            SSE_Solver, // parent class
            run_parallel, // func name
            
            );
        }

        int run_impl() override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            run_impl, // func name
            
            );
        }

        int init(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            init, // func name
            flag
            );
        }

        int final(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            final, // func name
            flag
            );
        }

        int cache(int flag) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Solver, // parent class
            cache, // func name
            flag
            );
        }
    };

    py::class_<SSE_Solver, Solver, PyTrampoline_SSE_Solver>(solvers_m, "SSE_Solver", py::dynamic_attr())
    .def(py::init<Param, Model*>())
    .def(py::init<const std::string&, Model*>())
    .def("ref_rand1", &SSE_Solver::ref_rand1, py::return_value_policy::reference_internal)
    .def("ref_rand2", &SSE_Solver::ref_rand2, py::return_value_policy::reference_internal)
    .def("ref_Eele", &SSE_Solver::ref_Eele, py::return_value_policy::reference_internal)
    .def("ref_DE", &SSE_Solver::ref_DE, py::return_value_policy::reference_internal)
    .def("ref_w", &SSE_Solver::ref_w, py::return_value_policy::reference_internal)
    .def("ref_tanhqwb", &SSE_Solver::ref_tanhqwb, py::return_value_policy::reference_internal)
    .def("ref_xicoeff1", &SSE_Solver::ref_xicoeff1, py::return_value_policy::reference_internal)
    .def("ref_xicoeff2", &SSE_Solver::ref_xicoeff2, py::return_value_policy::reference_internal)
    .def("ref_alpha_pref", &SSE_Solver::ref_alpha_pref, py::return_value_policy::reference_internal)
    .def("ref_CL", &SSE_Solver::ref_CL, py::return_value_policy::reference_internal)
    .def("ref_DEpW", &SSE_Solver::ref_DEpW, py::return_value_policy::reference_internal)
    .def("ref_DEmW", &SSE_Solver::ref_DEmW, py::return_value_policy::reference_internal)
    .def("ref_coeff_DEpW", &SSE_Solver::ref_coeff_DEpW, py::return_value_policy::reference_internal)
    .def("ref_coeff_DEmW", &SSE_Solver::ref_coeff_DEmW, py::return_value_policy::reference_internal)
    .def("ref_Hele", &SSE_Solver::ref_Hele, py::return_value_policy::reference_internal)
    .def("ref_Tele", &SSE_Solver::ref_Tele, py::return_value_policy::reference_internal)
    .def("ref_T0", &SSE_Solver::ref_T0, py::return_value_policy::reference_internal)
    .def("ref_Xzero", &SSE_Solver::ref_Xzero, py::return_value_policy::reference_internal)
    .def("ref_Xtran", &SSE_Solver::ref_Xtran, py::return_value_policy::reference_internal)
    .def("ref_Qzero", &SSE_Solver::ref_Qzero, py::return_value_policy::reference_internal)
    .def("ref_Qtran", &SSE_Solver::ref_Qtran, py::return_value_policy::reference_internal)
    .def("ref_BL", &SSE_Solver::ref_BL, py::return_value_policy::reference_internal)
    .def("ref_tmpbarr1", &SSE_Solver::ref_tmpbarr1, py::return_value_policy::reference_internal)
    .def("ref_tmpbarr2", &SSE_Solver::ref_tmpbarr2, py::return_value_policy::reference_internal)
    .def("ref_eac", &SSE_Solver::ref_eac, py::return_value_policy::reference_internal)
    .def("ref_eac_adia", &SSE_Solver::ref_eac_adia, py::return_value_policy::reference_internal)
    .def("ref_eac0", &SSE_Solver::ref_eac0, py::return_value_policy::reference_internal)
    .def("ref_rho0", &SSE_Solver::ref_rho0, py::return_value_policy::reference_internal)
    .def("ref_rhot", &SSE_Solver::ref_rhot, py::return_value_policy::reference_internal)
    .def("ref_fact1", &SSE_Solver::ref_fact1, py::return_value_policy::reference_internal)
    .def("ref_fact2", &SSE_Solver::ref_fact2, py::return_value_policy::reference_internal)
    .def("ref_fact3", &SSE_Solver::ref_fact3, py::return_value_policy::reference_internal)
    .def("ref_fact4", &SSE_Solver::ref_fact4, py::return_value_policy::reference_internal)
    .def("ref_eactmp", &SSE_Solver::ref_eactmp, py::return_value_policy::reference_internal)
    .def("ref_eactmp1", &SSE_Solver::ref_eactmp1, py::return_value_policy::reference_internal)
    .def("ref_eactmp2", &SSE_Solver::ref_eactmp2, py::return_value_policy::reference_internal)
    .def("ref_eacadd1", &SSE_Solver::ref_eacadd1, py::return_value_policy::reference_internal)
    .def("ref_eacadd2", &SSE_Solver::ref_eacadd2, py::return_value_policy::reference_internal)
    .def("ref_eactran", &SSE_Solver::ref_eactran, py::return_value_policy::reference_internal)
    .def("ref_crand1", &SSE_Solver::ref_crand1, py::return_value_policy::reference_internal)
    .def("ref_crand2", &SSE_Solver::ref_crand2, py::return_value_policy::reference_internal)
    .def("ref_invexpiEt", &SSE_Solver::ref_invexpiEt, py::return_value_policy::reference_internal)
    .def("ref_U", &SSE_Solver::ref_U, py::return_value_policy::reference_internal)
    .def("ref_expiwt", &SSE_Solver::ref_expiwt, py::return_value_policy::reference_internal)
    .def("ref_expiDEt", &SSE_Solver::ref_expiDEt, py::return_value_policy::reference_internal)
    .def("ref_invexpiwt", &SSE_Solver::ref_invexpiwt, py::return_value_policy::reference_internal)
    .def("ref_invexpiDEt", &SSE_Solver::ref_invexpiDEt, py::return_value_policy::reference_internal)
    .def("ref_expiwdht", &SSE_Solver::ref_expiwdht, py::return_value_policy::reference_internal)
    .def("ref_invexpiwdht", &SSE_Solver::ref_invexpiwdht, py::return_value_policy::reference_internal)
    .def("ref_expiwt_now", &SSE_Solver::ref_expiwt_now, py::return_value_policy::reference_internal)
    .def("ref_invexpiwt_now", &SSE_Solver::ref_invexpiwt_now, py::return_value_policy::reference_internal)
    .def("ref_Xtmpj", &SSE_Solver::ref_Xtmpj, py::return_value_policy::reference_internal)
    .def("ref_Xtmp2j", &SSE_Solver::ref_Xtmp2j, py::return_value_policy::reference_internal)
    .def("ref_xi", &SSE_Solver::ref_xi, py::return_value_policy::reference_internal)
    .def("ref_time_ker", &SSE_Solver::ref_time_ker, py::return_value_policy::reference_internal)
    .def("ref_mem_arr", &SSE_Solver::ref_mem_arr, py::return_value_policy::reference_internal)
    .def_static("name", &SSE_Solver::name)
    .def("init", static_cast<int (SSE_Solver::*)(const int&)>(&SSE_Solver::init))
    .def("correlation", static_cast<int (SSE_Solver::*)(const int&, NAD_TCFer&)>(&SSE_Solver::correlation))
    .def("sampler", static_cast<int (SSE_Solver::*)(const int&, NAD_TCFer&)>(&SSE_Solver::sampler))
    .def("call_memory_array_adia", static_cast<int (SSE_Solver::*)(const num_real&)>(&SSE_Solver::call_memory_array_adia))
    .def("sse", static_cast<int (SSE_Solver::*)(NAD_TCFer&, const int&, const int&)>(&SSE_Solver::sse))
    .def("sse_adia", static_cast<int (SSE_Solver::*)(NAD_TCFer&, const int&, const int&)>(&SSE_Solver::sse_adia))
    .def("run_parallel", static_cast<int (SSE_Solver::*)()>(&SSE_Solver::run_parallel))
    .def("run", static_cast<int (Solver::*)()>(&Solver::run))
    .def("run_impl", static_cast<int (Solver::*)()>(&Solver::run_impl))
    .def("run_parallel", static_cast<int (Solver::*)()>(&Solver::run_parallel))
    .def("init", static_cast<int (Solver::*)(int)>(&Solver::init))
    .def("final", static_cast<int (Solver::*)(int)>(&Solver::final))
    .def("cache", static_cast<int (Solver::*)(int)>(&Solver::cache));

    class PyTrampoline_Bath : public Bath {
        public:
        using Bath::Bath;
    };

    py::class_<Bath, Model, PyTrampoline_Bath>(models_m, "Bath", py::dynamic_attr())
    .def(py::init<const Param&, const int&, const int&>())
    .def(py::init<const std::string&, const int&, const int&>())
    .def_readwrite("coup_type", &Bath::coup_type)
    .def_readwrite("coup_flag", &Bath::coup_flag)
    .def_readwrite("spec_type", &Bath::spec_type)
    .def_readwrite("spec_flag", &Bath::spec_flag)
    .def_readwrite("omegac", &Bath::omegac)
    .def_readwrite("lambda", &Bath::lambda)
    .def_readwrite("beta", &Bath::beta)
    .def_readwrite("nbath", &Bath::nbath)
    .def_readwrite("F", &Bath::F)
    .def("ref_Q", &Bath::ref_Q, py::return_value_policy::reference_internal)
    .def("Discretization", [](Bath& self, py::array_t<double, py::array::c_style | py::array::forcecast> W_arr_arr, 
        py::array_t<double, py::array::c_style | py::array::forcecast> C_arr_arr, 
        const int& Nb) {
            return self.Discretization(W_arr_arr.mutable_data(), C_arr_arr.mutable_data(), Nb); 
        }
    )
    .def("J_Ohmic", static_cast<double (Bath::*)(const double&)>(&Bath::J_Ohmic))
    .def("J_Debye", static_cast<double (Bath::*)(const double&)>(&Bath::J_Debye))
    .def("J_Refit", [](Bath& self, const double& w, 
        py::array_t<double, py::array::c_style | py::array::forcecast> W_arr_arr, 
        py::array_t<double, py::array::c_style | py::array::forcecast> C_arr_arr, 
        const int& Nb) {
            return self.J_Refit(w, W_arr_arr.mutable_data(), C_arr_arr.mutable_data(), Nb); 
        }
    )
    .def("J", [](Bath& self, const double& w, 
        py::array_t<double, py::array::c_style | py::array::forcecast> W_arr_arr, 
        py::array_t<double, py::array::c_style | py::array::forcecast> C_arr_arr, 
        const int& Nb) {
            return self.J(w, W_arr_arr.mutable_data(), C_arr_arr.mutable_data(), Nb); 
        }
    )
    .def("fun_Cw_Re", [](Bath& self, const double& w, 
        py::array_t<double, py::array::c_style | py::array::forcecast> W_arr_arr, 
        py::array_t<double, py::array::c_style | py::array::forcecast> C_arr_arr, 
        const int& Nb) {
            return self.fun_Cw_Re(w, W_arr_arr.mutable_data(), C_arr_arr.mutable_data(), Nb); 
        }
    );

    class PyTrampoline_Thermostat : public Thermostat {
        public:
        using Thermostat::Thermostat;

        int evolve(num_real* nr, num_real* np, num_real* nm, const num_real& dt, const int& N, const int& start = 0,
                       num_real gammal_in = -1.0f) override {
            PYBIND11_OVERRIDE(
            int, // return type
            Thermostat, // parent class
            evolve, // func name
            nr, np, nm, dt, N, start, gammal_in
            );
        }
    };

    py::class_<Thermostat, Model, PyTrampoline_Thermostat>(models_m, "Thermostat", py::dynamic_attr())
    .def(py::init<const Param&>())
    .def(py::init<const std::string&>())
    .def_readwrite("beta", &Thermostat::beta)
    .def("ref_nhc_r", &Thermostat::ref_nhc_r, py::return_value_policy::reference_internal)
    .def("ref_nhc_p", &Thermostat::ref_nhc_p, py::return_value_policy::reference_internal)
    .def("ref_nhc_G", &Thermostat::ref_nhc_G, py::return_value_policy::reference_internal)
    .def("ref_nhc_Q", &Thermostat::ref_nhc_Q, py::return_value_policy::reference_internal)
    .def("init_alloc", static_cast<int (Thermostat::*)(const int&)>(&Thermostat::init_alloc))
    .def("dothermo", static_cast<bool (Thermostat::*)(const int&)>(&Thermostat::dothermo))
    .def("set_gammal", static_cast<void (Thermostat::*)(num_real)>(&Thermostat::set_gammal))
    .def("evolve", [](Thermostat& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const num_real& dt, 
        const int& N, 
        const int& start, 
        num_real gammal_in) {
            return self.evolve(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), dt, N, start, gammal_in); 
        }
    )
    .def("run_NHC", [](Thermostat& self, py::array_t<num_real, py::array::c_style | py::array::forcecast> nr_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> np_arr, 
        py::array_t<num_real, py::array::c_style | py::array::forcecast> nm_arr, 
        const num_real& dt, 
        const int& N, 
        const int& start, 
        num_real gammal_in) {
            return self.run_NHC(nr_arr.mutable_data(), np_arr.mutable_data(), nm_arr.mutable_data(), dt, N, start, gammal_in); 
        }
    );
}
// clang-format off


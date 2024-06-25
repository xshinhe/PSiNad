py::class_<Param, std::shared_ptr<Param>> PyParam(m, "Param", py::dynamic_attr());

py::enum_<Param::LoadOption>(PyParam, "LoadOption")
    .value("fromString", Param::LoadOption::fromString)
    .value("fromFile", Param::LoadOption::fromFile)
    .export_values();

PyParam.def(py::init<const std::string&, Param::LoadOption>());

PyParam.def("has_key", &Param::has_key, py::return_value_policy::reference_internal);

PyParam.def(
    "get_bool", [](Param& self, const std::vector<std::string>& keys) { return self.get_bool(keys); },
    py::return_value_policy::reference_internal);

PyParam.def(
    "get_int", [](Param& self, const std::vector<std::string>& keys) { return self.get_int(keys); },
    py::return_value_policy::reference_internal);

PyParam.def(
    "get_real",
    [](Param& self, const std::vector<std::string>& keys, phys::dimension7 qdim) {
        return self.get_double(keys, "__loc__", qdim);
    },
    py::arg("keys"),                 //
    py::arg("qdim") = phys::none_d,  //
    py::return_value_policy::reference_internal);

PyParam.def(
    "get_str", [](Param& self, const std::vector<std::string>& keys) { return self.get_string(keys); },
    py::return_value_policy::reference_internal);


// PyParam.def("get_int", &Param::get<int,true>, py::arg("key") = "", py::arg("loc") = "",
//    py::arg("qdim") = phys::none_d,
//    py::arg("default") = int());
// PyParam.def("get_real", &Param::get<double,true>, py::arg("key") = "", py::arg("loc") = "",
//    py::arg("qdim") = phys::none_d,
//    py::arg("default") = double());
// PyParam.def("get_bool", &Param::get<bool,true>, py::arg("key") = "", py::arg("loc") = "",
//    py::arg("qdim") = phys::none_d,
//    py::arg("default") = bool());

// PyParam.def("get_str", &Param::get<std::string,true>, py::arg("key") = "", py::arg("loc") = "",
//     py::arg("qdim") = phys::none_d,
//     py::arg("default") = "");

PyParam.def("__repr__", &Param::repr);

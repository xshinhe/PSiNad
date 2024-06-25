py::class_<Param, std::shared_ptr<RuleSet>>(m, "RuleSet", py::dynamic_attr())
    .def(py::init<const std::string &, Param::LoadOption>())
    .def("has_key", &Param::has_key, py::return_value_policy::reference_internal);

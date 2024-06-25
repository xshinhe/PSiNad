py::class_<Status, std::shared_ptr<Status>>(m, "Status", py::dynamic_attr())
    .def(py::init<>())
    .def(py::init<bool, int, int, int>())
    .def(py::init<const Status &>());
// .def("__repr__", [](Status& self){});

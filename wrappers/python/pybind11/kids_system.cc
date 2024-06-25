py::class_<System, std::shared_ptr<System>>(m, "System")
    .def(py::init<std::shared_ptr<Model>, std::shared_ptr<Param>, std::shared_ptr<DataSet>>());

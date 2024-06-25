py::class_<Solver, std::shared_ptr<Solver>>(m, "Solver", py::dynamic_attr())  //
    .def(py::init<std::shared_ptr<Kernel>>())
    .def("getSystem", &Solver::getSystem)
    .def("getSolverKernel", &Solver::getSolverKernel);

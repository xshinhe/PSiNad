py::class_<Solver, std::shared_ptr<Solver>> PySolver(m, "Solver", py::dynamic_attr());

PySolver.def(py::init<std::shared_ptr<Kernel>>());

PySolver.def("getSystem", &Solver::getSystem)
    .def("getSolverKernel", &Solver::getSolverKernel)
    .def("setInputParam", &Solver::setInputParam)
    .def("setInputDataSet", &Solver::setInputDataSet);

PySolver.def("addApplication", [](Solver& self, Kernel& K) {
        int a = 0;
        return self;
    });

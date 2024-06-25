class PyKernel : public Kernel {
   public:
    /* Inherit the constructors */
    using Kernel::Kernel;

    /* Trampoline (need one for each virtual function) */
    const std::string getName() override {
        PYBIND11_OVERRIDE_PURE(const std::string, /* Return type */
                               Kernel,            /* Parent class */
                               getName,           /* Name of function in C++ (must match Python name) */
        );
    }

    void setInputParam_impl(std::shared_ptr<Param> PM) override {
        PYBIND11_OVERRIDE_PURE(void,               /* Return type */
                               Kernel,             /* Parent class */
                               setInputParam_impl, /* Name of function in C++ (must match Python name) */
                               PM);
    }

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS) override {
        PYBIND11_OVERRIDE_PURE(void,                 /* Return type */
                               Kernel,               /* Parent class */
                               setInputDataSet_impl, /* Name of function in C++ (must match Python name) */
                               DS);
    }

    Status& initializeKernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(Status&,               /* Return type */
                               Kernel,                /* Parent class */
                               initializeKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
    Status& executeKernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(Status&,            /* Return type */
                               Kernel,             /* Parent class */
                               executeKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
};

py::class_<Kernel, PyKernel, std::shared_ptr<Kernel>>(m, "Kernel", py::dynamic_attr())  //
    .def(py::init<const std::string&>())
    .def("getName", &Kernel::getName)
    .def("appendChild", &Kernel::appendChild)
    .def("insertAt", &Kernel::insertAt)
    .def("removeAt", &Kernel::removeAt)
    .def("updateAt", &Kernel::updateAt)
    .def("generateInformationString", &Kernel::generateInformationString)
    .def("setTiming", &Kernel::setTiming)
    .def("setInputParam", &Kernel::setInputParam)
    .def("setInputDataSet", &Kernel::setInputDataSet)
    .def("setRuleSet", &Kernel::setRuleSet)
    .def("initializeKernel", &Kernel::initializeKernel)
    .def("executeKernel", &Kernel::executeKernel)
    .def("finalizeKernel", &Kernel::finalizeKernel);

m.def("defaultModelFactory", &defaultModelFactory);


// template <class T>
// defaultSolverFactory_T(const std::string&, T);

m.def("defaultSolverFactory",
      static_cast<std::shared_ptr<Solver> (*)(const std::string&, std::shared_ptr<Model>)>(&defaultSolverFactory));
m.def("defaultSolverFactory",
      static_cast<std::shared_ptr<Solver> (*)(const std::string&, std::shared_ptr<System>)>(&defaultSolverFactory));

// m.def("defaultSolverFactory", py::overload_cast<const std::string&, std::shared_ptr<System>>(&defaultSolverFactory));

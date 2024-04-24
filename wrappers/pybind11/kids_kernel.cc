class PyKernel : public Kernel {
   public:
    /* Inherit the constructors */
    using Kernel::Kernel;

    /* Trampoline (need one for each virtual function) */
    const std::string name() override {
        PYBIND11_OVERRIDE_PURE(const std::string, /* Return type */
                               Kernel,            /* Parent class */
                               name,              /* Name of function in C++ (must match Python name) */
        );
    }

    void inputParam_impl(Param* PM) override {
        PYBIND11_OVERRIDE_PURE(void,            /* Return type */
                               Kernel,          /* Parent class */
                               inputParam_impl, /* Name of function in C++ (must match Python name) */
                               PM);
    }

    void inputDataSet_impl(DataSet* DS) override {
        PYBIND11_OVERRIDE_PURE(void,              /* Return type */
                               Kernel,            /* Parent class */
                               inputDataSet_impl, /* Name of function in C++ (must match Python name) */
                               DS);
    }

    Status& initKernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(void,            /* Return type */
                               Kernel,          /* Parent class */
                               initKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
    Status& execkernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(int,             /* Return type */
                               Kernel,          /* Parent class */
                               execKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
};

py::class_<Kernel, PyKernel, std::shared_ptr<Kernel>>(m, "Kernel", py::dynamic_attr())  //
    .def(py::init<const std::string&>())
    .def("name", &Kernel::name)
    .def("push", &Kernel::push)
    .def("insert", &Kernel::insert)
    .def("erase", &Kernel::erase)
    .def("update", &Kernel::update)
    .def("scheme", &Kernel::scheme)
    .def("inputParam", &Kernel::inputParam)
    .def("inputDataSet", &Kernel::inputDataSet)
    .def("initKernel", &Kernel::initKernel)
    .def("execKernel", &Kernel::execKernel)
    .def("finalKernel", &Kernel::finalKernel);

m.def("modelfactory", &ModelFactory);

m.def("solverfactory", &SolverFactory);

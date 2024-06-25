class PyModel : public Model {
   public:
    /* Inherit the constructors */
    using Model::Model;

    /* Trampoline (need one for each virtual function) */
    const std::string getName() override {
        PYBIND11_OVERRIDE_PURE(const std::string, /* Return type */
                               Model,             /* Parent class */
                               getName,           /* Name of function in C++ (must match Python name) */
        );
    }

    void setInputParam_impl(std::shared_ptr<Param> PM) override {
        PYBIND11_OVERRIDE_PURE(void,               /* Return type */
                               Model,              /* Parent class */
                               setInputParam_impl, /* Name of function in C++ (must match Python name) */
                               PM);
    }

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS) override {
        PYBIND11_OVERRIDE_PURE(void,                 /* Return type */
                               Model,                /* Parent class */
                               setInputDataSet_impl, /* Name of function in C++ (must match Python name) */
                               DS);
    }

    Status& initializeKernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(Status&,               /* Return type */
                               Model,                 /* Parent class */
                               initializeKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
    Status& executeKernel_impl(Status& stat) override {
        PYBIND11_OVERRIDE_PURE(Status&,            /* Return type */
                               Model,              /* Parent class */
                               executeKernel_impl, /* Name of function in C++ (must match Python name) */
                               stat);
    }
};

py::class_<Model, PyModel, std::shared_ptr<Model>>(m, "Model", py::dynamic_attr())  //
    .def(py::init<const std::string&>());

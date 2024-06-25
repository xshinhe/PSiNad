py::class_<Platform, std::shared_ptr<Platform>>(m, "Platform", py::dynamic_attr())  //
    .def("getName", &Platform::getName);                                            //

m.def("usePlatform", &Platform::usePlatform);

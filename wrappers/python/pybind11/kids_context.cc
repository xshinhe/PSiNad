py::class_<Context> PyContext(m, "Context", py::dynamic_attr());

PyContext.def(py::init<std::shared_ptr<Platform>, std::shared_ptr<System>, std::vector<std::shared_ptr<Solver>>>());

PyContext.def("run", [](Context& self, Status& stat) {
    self.run(stat);
    return stat;
});

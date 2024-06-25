py::class_<Context> PyContext(m, "Context", py::dynamic_attr());

PyContext.def(py::init<std::shared_ptr<Platform>, std::shared_ptr<System>, std::vector<std::shared_ptr<Solver>>>());

PyContext.def("execute", [](Context& self) {
    Status stat;
    self.execute(stat);
    return stat.succ;
});

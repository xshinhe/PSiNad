py::class_<Shape> PyShape(m, "Shape", py::dynamic_attr());

PyShape.def(py::init<std::vector<std::size_t>>()).def(py::init<std::size_t>());

py::class_<DataSet> PyDataSet(m, "DataSet", py::dynamic_attr());

PyDataSet.def(py::init<>());

PyDataSet.def(
    "_def",
    [](DataSet& self, const std::string& key, py::tuple pyTuple, const std::string& dtype, const std::string& doc) {
        std::vector<std::size_t> shape;
        shape.reserve(pyTuple.size());  // Reserve space for efficiency
        for (size_t i = 0; i < pyTuple.size(); ++i) {
            int value = pyTuple[i].cast<int>();
            shape.push_back(value);
        }
        if (dtype == "int") {
            self._def_int(key, Shape{shape}, doc);
        } else if (dtype == "real") {
            self._def_real(key, Shape{shape}, doc);
        } else if (dtype == "complex") {
            self._def_complex(key, Shape{shape}, doc);
        } else {
            // throw std::runtime_error("can not define with this dtype!");
        }
    },
    py::arg("key"), py::arg("shape"), py::arg("dtype") = "int", py::arg("doc") = "");

PyDataSet.def("_undef", &DataSet::_undef, py::arg("key") = "0");

PyDataSet.def("numpy", [](DataSet& self, const std::string& key) {
    auto       inode  = self.node(key);
    kids_dtype n_type = inode->type();
    switch (n_type) {
        case kids_int_type: {
            std::size_t n_size = static_cast<Tensor<kids_int>*>(inode)->size();
            void*       n_data = static_cast<Tensor<kids_int>*>(inode)->data();
            return py::array({n_size},                                          // shape
                             {sizeof(int)},                                     // stride
                             (int*) n_data,                                     // data pointer
                             py::capsule(n_data, [](void* _void_n_data) { ; })  // zero-copy cost
            );
            break;
        }
        case kids_real_type: {
            std::size_t n_size = static_cast<Tensor<kids_real>*>(inode)->size();
            void*       n_data = static_cast<Tensor<kids_real>*>(inode)->data();
            return py::array({n_size},                                          // shape
                             {sizeof(kids_real)},                               // stride
                             (kids_real*) n_data,                               // data pointer
                             py::capsule(n_data, [](void* _void_n_data) { ; })  // zero-copy cost
            );
            break;
        }
        case kids_complex_type: {
            std::size_t n_size = static_cast<Tensor<kids_complex>*>(inode)->size();
            void*       n_data = static_cast<Tensor<kids_complex>*>(inode)->data();
            return py::array({n_size},                                          // shape
                             {sizeof(kids_complex)},                            // stride
                             (kids_complex*) n_data,                            // data pointer
                             py::capsule(n_data, [](void* _void_n_data) { ; })  // zero-copy cost
            );
            break;
        }
        default:
            throw std::runtime_error("can not converted to numpy!");
    }
    return py::array();
});

PyDataSet.def(
    "help", [](DataSet& self, const std::string& name) { return self.help(name); }, py::arg("name") = "");

PyDataSet.def("__repr__", [](DataSet& self) { return self.repr(); });

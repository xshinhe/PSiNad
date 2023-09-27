py::buffer_info(m.data(),                               /* Pointer to buffer */
                sizeof(float),                          /* Size of one scalar */
                py::format_descriptor<float>::format(), /* Python struct-style format descriptor */
                2,                                      /* Number of dimensions */
                {m.rows(), m.cols()},                   /* Buffer dimensions */
                {sizeof(float) * m.cols(),              /* Strides (in bytes) for each index */
                 sizeof(float)});

struct buffer_info {
    void* ptr;
    py::ssize_t itemsize;
    std::string format;
    py::ssize_t ndim;
    std::vector<py::ssize_t> shape;
    std::vector<py::ssize_t> strides;
};

enum StateType {
    StateNone,
    StateBool,
    StateInt,
    StateReal,
    StateComplex,
};

template <typename T>
class StateHelper {
   public:
    StateType type = StateNone;
};
template <>
class StateHelper<bool> {
   public:
    StateType type = StateBool;
};
template <>
class StateHelper<int> {
   public:
    StateType type = StateInt;
};
template <>
class StateHelper<double> {
   public:
    StateType type = StateReal;
};
template <>
class StateHelper<std::complex<double>> {
   public:
    StateType type = StateComplex;
};



class State {
   public:
    State(double* ptr, std::size_t size);
    template <typename T>
    State(std::vector<T> vec);

   protected:
    std::string name;
    int size;

    vector<int> _int;
    vector<double> _double;
    vector<std::complex<double>> _complex;
};
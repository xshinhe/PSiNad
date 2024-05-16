#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kids/DataSet.h"
#include "kids/Kernel.h"
#include "kids/ModelFactory.h"
#include "kids/Param.h"
#include "kids/SolverFactory.h"

namespace py = pybind11;
using namespace PROJECT_NS;

PYBIND11_MODULE(libpykids_v1, m) {
    // clang-format off
#include "kids_phys.cc"
#include "kids_param.cc"
#include "kids_dataset.cc"
#include "kids_kernel.cc"
    // clang-format on
}

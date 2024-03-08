#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/DataSet.h"
#include "../core/Kernel.h"
#include "../core/Param.h"
#include "../models/ModelFactory.h"
#include "../solvers/SolverFactory.h"

namespace py = pybind11;
using namespace PROJECT_NS;

PYBIND11_MODULE(libpykids, m) {
#include "kids_phys.bind"
//
#include "kids_param.bind"
//
#include "kids_dataset.bind"
//
#include "kids_kernel.bind"
}

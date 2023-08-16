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

PYBIND11_MODULE(libopendf, m) {
#include "opendf_phys.bind"
//
#include "opendf_param.bind"
//
#include "opendf_dataset.bind"
//
#include "opendf_kernel.bind"
}

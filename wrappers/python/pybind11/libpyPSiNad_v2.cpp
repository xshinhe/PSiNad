#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "kids/Context.h"
#include "kids/DataSet.h"
#include "kids/Kernel.h"
#include "kids/Model.h"
#include "kids/ModelFactory.h"
#include "kids/Param.h"
#include "kids/Platform.h"
#include "kids/Solver.h"
#include "kids/SolverFactory.h"
#include "kids/System.h"
#include "kids/Types.h"
#include "kids/Variable.h"
#include "kids/chem.h"
#include "kids/phys.h"
#include "kids/vars_list.h"

namespace py = pybind11;
using namespace PROJECT_NS;

PYBIND11_MODULE(libpyPSiNad_v2, m) {
    // clang-format off
#include "kids_phys.cc"
#include "kids_chem.cc"
#include "kids_status.cc"
#include "kids_var.cc"
#include "kids_param.cc"
#include "kids_dataset.cc"
#include "kids_kernel.cc"
#include "kids_model.cc"
#include "kids_system.cc"
#include "kids_solver.cc"
#include "kids_context.cc"
#include "kids_platform.cc"
    // clang-format on
}

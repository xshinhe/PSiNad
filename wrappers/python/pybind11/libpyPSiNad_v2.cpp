#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "psnd/Context.h"
#include "psnd/DataSet.h"
#include "psnd/Kernel.h"
#include "psnd/Model.h"
#include "psnd/ModelFactory.h"
#include "psnd/Param.h"
#include "psnd/Platform.h"
#include "psnd/Solver.h"
#include "psnd/SolverFactory.h"
#include "psnd/System.h"
#include "psnd/Types.h"
#include "psnd/Variable.h"
#include "psnd/chem.h"
#include "psnd/phys.h"
#include "psnd/vars_list.h"

namespace py = pybind11;
using namespace PROJECT_NS;

PYBIND11_MODULE(libpyPSiNaD_v2, m) {
    // clang-format off
#include "psnd_phys.cc"
#include "psnd_chem.cc"
#include "psnd_status.cc"
#include "psnd_var.cc"
#include "psnd_param.cc"
#include "psnd_dataset.cc"
#include "psnd_kernel.cc"
#include "psnd_model.cc"
#include "psnd_system.cc"
#include "psnd_solver.cc"
#include "psnd_context.cc"
#include "psnd_platform.cc"
    // clang-format on
}

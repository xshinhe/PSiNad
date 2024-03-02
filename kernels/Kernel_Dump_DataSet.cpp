#include "Kernel_Dump_DataSet.h"

#include "../core/linalg.h"

namespace kids {

void Kernel_Dump_DataSet::read_param_impl(Param* PM) { fn = PM->get<std::string>("dump", LOC(), "final"); }

int Kernel_Dump_DataSet::exec_kernel_impl(int stat) {
    if (fn == "null") return 0;
    try {
        _DataSet->dump(utils::concat(fn, stat, ".ds"));
    } catch (std::runtime_error& e) { throw basic_error(fn); }
    return 0;
}

};  // namespace kids

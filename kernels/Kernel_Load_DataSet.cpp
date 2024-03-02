#include "Kernel_Load_DataSet.h"

#include "../core/linalg.h"

namespace kids {

void Kernel_Load_DataSet::read_param_impl(Param* PM) { fn = PM->get<std::string>("load", LOC(), "NULL"); }

int Kernel_Load_DataSet::exec_kernel_impl(int stat) {
    if (fn != "NULL") {
        try {
            _DataSet->load(fn);
        } catch (std::runtime_error& e) { throw basic_error(fn); }
    }
    return 0;
}

};  // namespace kids

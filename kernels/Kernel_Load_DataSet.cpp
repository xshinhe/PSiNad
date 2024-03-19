#include "Kernel_Load_DataSet.h"

#include "../core/linalg.h"

namespace PROJECT_NS {

void Kernel_Load_DataSet::read_param_impl(Param* PM) { fn = PM->get<std::string>("load", LOC(), "NULL"); }

int Kernel_Load_DataSet::exec_kernel_impl(int stat) {
    if (fn == "" || fn == "NULL" || fn == "null") return 0;
    try {
        std::ifstream ifs{fn};
        _DataSet->load(ifs);
        ifs.close();
    } catch (std::runtime_error& e) { throw basic_error(fn); }
    return 0;
}

};  // namespace PROJECT_NS

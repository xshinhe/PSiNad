#include "kids/Kernel_Load_DataSet.h"

#include "kids/linalg.h"

namespace PROJECT_NS {

void Kernel_Load_DataSet::setInputParam_impl(std::shared_ptr<Param>& PM) { fn = PM->get_string("load", LOC(), "NULL"); }

Status& Kernel_Load_DataSet::executeKernel_impl(Status& stat) {
    if (fn == "" || fn == "NULL" || fn == "null") return 0;
    try {
        std::ifstream ifs{fn};
        _dataset->load(ifs);
        ifs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return 0;
}

};  // namespace PROJECT_NS

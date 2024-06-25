#include "kids/Kernel_Load_DataSet.h"

#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

const std::string Kernel_Load_DataSet::getName() { return "Kernel_Load_DataSet"; }

int Kernel_Load_DataSet::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Load_DataSet::setInputParam_impl(std::shared_ptr<Param> PM) { fn = PM->get_string("load", LOC(), "NULL"); }

Status& Kernel_Load_DataSet::executeKernel_impl(Status& stat) {
    if (fn == "" || fn == "NULL" || fn == "null") return stat;
    try {
        std::ifstream ifs{fn};
        _dataset->load(ifs);
        ifs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return stat;
}

};  // namespace PROJECT_NS

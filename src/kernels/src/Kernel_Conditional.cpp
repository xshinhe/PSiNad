#include "kids/Kernel_Conditional.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Conditional::getName() { return "Kernel_Conditional"; }

int Kernel_Conditional::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Conditional::setInputParam_impl(std::shared_ptr<Param> PM){};

void Kernel_Conditional::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    at_condition = DS->def(DATA::flowcontrol::at_condition);
}

Status& Kernel_Conditional::initializeKernel_impl(Status& stat) {
    at_condition[0] = false;
    return stat;
}

Status& Kernel_Conditional::executeKernel_impl(Status& stat) {
    if (at_condition[0]) {
        for (auto& pkernel : _child_kernels) { pkernel->executeKernel(stat); }
    }
    return stat;
}
};  // namespace PROJECT_NS

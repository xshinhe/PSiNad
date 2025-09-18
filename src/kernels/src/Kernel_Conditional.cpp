#include "psnd/Kernel_Conditional.h"

#include "psnd/hash_fnv1a.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Conditional::getName() { return "Kernel_Conditional"; }

int Kernel_Conditional::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Conditional::setInputParam_impl(std::shared_ptr<Param> PM) {};

void Kernel_Conditional::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    at_condition = DS->def(DATA::control::at_condition);
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

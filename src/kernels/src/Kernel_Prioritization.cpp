#include "kids/Kernel_Prioritization.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel_Prioritization::Kernel_Prioritization(std::vector<std::shared_ptr<Kernel>> kers, int ptype_in)
    : Kernel(), ptype{ptype_in} {
    for (auto& ker : kers) { _ref_kernels.push_back(ker); }
}

const std::string Kernel_Prioritization::getName() {
    std::stringstream ss;
    ss << "Kernel_Prioritization [" << ptype << "]";
    for (auto& ker : _ref_kernels) ss << " #" << std::setfill('0') << std::setw(2) << ker->getID();
    return ss.str();
}

int Kernel_Prioritization::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Prioritization::setInputParam_impl(std::shared_ptr<Param>& PM) {
    if (ptype == 0)
        for (auto& ker : _ref_kernels) ker->setInputParam(PM);
}

void Kernel_Prioritization::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    if (ptype == 1)
        for (auto& ker : _ref_kernels) ker->setInputDataSet(DS);
}

Status& Kernel_Prioritization::initializeKernel_impl(Status& stat) {
    if (ptype == 2)
        for (auto& ker : _ref_kernels) ker->initializeKernel(stat);
    return stat;
}

};  // namespace PROJECT_NS
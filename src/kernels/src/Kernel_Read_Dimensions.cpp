#include "kids/Kernel_Read_Dimensions.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel_Read_Dimensions::Kernel_Read_Dimensions() : Kernel(){};

const std::string Kernel_Read_Dimensions::getName() { return "Kernel_Read_Dimensions"; }

int Kernel_Read_Dimensions::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Read_Dimensions::setInputParam_impl(std::shared_ptr<Param> PM) {
    Dimension::M = _param->get_int({"solver.M"}, LOC(), 1);
    Dimension::P = _param->get_int({"solver.P"}, LOC(), 1);
    Dimension::N = _param->get_int({"model.N"}, LOC(), 1);
    Dimension::F = _param->get_int({"model.F"}, LOC(), 1);
    Dimension::static_build_shapes();
};

};  // namespace PROJECT_NS

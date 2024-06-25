#include "kids/Kernel_Read_Dimensions.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel_Read_Dimensions::Kernel_Read_Dimensions() : Kernel(){};

const std::string Kernel_Read_Dimensions::getName() { return "Kernel_Read_Dimensions"; }

int Kernel_Read_Dimensions::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Read_Dimensions::setInputParam_impl(std::shared_ptr<Param> PM) {
    Dimension::M = PM->get_int("M", LOC(), 1);
    Dimension::P = PM->get_int("P", LOC(), 1);
    Dimension::N = PM->get_int("N", LOC(), 1);
    Dimension::F = PM->get_int("F", LOC(), 1);
    Dimension::static_build_shapes();
};

};  // namespace PROJECT_NS

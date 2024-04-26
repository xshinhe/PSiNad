#include "kids/Kernel_Custom.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel_Custom::Kernel_Custom(const std::string& customized_name) : Kernel(customized_name){};

const std::string Kernel_Custom::getName() { return kernel_name; }

int Kernel_Custom::getType() const { return utils::hash(FUNCTION_NAME); }

};  // namespace PROJECT_NS
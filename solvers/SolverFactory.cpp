#include "SolverFactory.h"

namespace PROJECT_NS {

extern std::shared_ptr<Kernel> NAD_Kernel(std::shared_ptr<Kernel> kmodel, std::string NAD_Kernel_name);

std::shared_ptr<Kernel> SolverFactory(const std::string& name, std::shared_ptr<Kernel> kmodel) {
    if (false) {
    } else if (name == "Hello") {
        // return Hello_SBuilder(kmodel);
    } else if (name == "CMM" || name == "SQC" || name == "MMD" || name == "SH" || name == "MMSH" || name == "MCE") {
        return NAD_Kernel(kmodel, name);
    } else {
        throw std::runtime_error("unknown solver name");
    }
    return nullptr;
}

};  // namespace PROJECT_NS

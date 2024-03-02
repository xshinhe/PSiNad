#include "SolverFactory.h"

namespace kids {

extern std::shared_ptr<Kernel> NAD_Kernel(std::shared_ptr<Kernel> kmodel, std::string NAD_Kernel_name);
extern std::shared_ptr<Kernel> NAD_Adapt_Kernel(std::shared_ptr<Kernel> kmodel, std::string NAD_Kernel_name);

std::shared_ptr<Kernel> SolverFactory(const std::string& name, std::shared_ptr<Kernel> kmodel) {
    if (false) {
    } else if (name == "Hello") {
        // return Hello_SBuilder(kmodel);
    } else if (name == "CMM" || name == "SQC" || name == "MMD" || name == "SH" || name == "CMSH" || name == "MMSH" ||
               name == "MCE") {
        // return NAD_Kernel(kmodel, name);
        return NAD_Adapt_Kernel(kmodel, name);
    } else {
        throw std::runtime_error("unknown solver name");
    }
    return nullptr;
}

};  // namespace kids

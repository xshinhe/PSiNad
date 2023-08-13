#include "SolverFactory.h"

namespace PROJECT_NS {

extern std::shared_ptr<Kernel> CMM_Kernel(std::shared_ptr<Kernel> kmodel);
extern std::shared_ptr<Kernel> SQC_Kernel(std::shared_ptr<Kernel> kmodel);
extern std::shared_ptr<Kernel> MMD_Kernel(std::shared_ptr<Kernel> kmodel);
extern std::shared_ptr<Kernel> SH_Kernel(std::shared_ptr<Kernel> kmodel);
extern std::shared_ptr<Kernel> MMSH_Kernel(std::shared_ptr<Kernel> kmodel);

std::shared_ptr<Kernel> SolverFactory(const std::string& name, std::shared_ptr<Kernel> kmodel) {
    if (false) {
    } else if (name == "Hello") {
        // return Hello_SBuilder(kmodel);
    } else if (name == "CMM") {
        return CMM_Kernel(kmodel);
    } else if (name == "SQC") {
        return SQC_Kernel(kmodel);
    } else if (name == "MMD") {
        return MMD_Kernel(kmodel);
    } else if (name == "SH") {
        return SH_Kernel(kmodel);
    } else if (name == "MMSH") {
        return MMSH_Kernel(kmodel);
    } else {
        throw std::runtime_error("unknown solver name");
    }
    return nullptr;
}

};  // namespace PROJECT_NS

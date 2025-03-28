#include "kids/SolverFactory.h"

namespace PROJECT_NS {

extern std::shared_ptr<Solver> Sampling_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name);

extern std::shared_ptr<Solver> NAD_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name);

extern std::shared_ptr<Solver> NAD_Coup_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name);

extern std::shared_ptr<Solver> NAD_Adapt_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name);

extern std::shared_ptr<Solver> NAD_AdaptM_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name);

std::shared_ptr<Solver> defaultSolverFactory(const std::string& name, std::shared_ptr<Model> kmodel) {
    if (false) {
    } else if (name == "Hello") {
        // return Hello_SBuilder(kmodel);
    } else if (name == "Sampling") {
        return Sampling_Kernel(kmodel, name);
    } else if (name == "NAD") {
        return NAD_Kernel(kmodel, name);
    } else if (name == "NAD-coup") {
        return NAD_Coup_Kernel(kmodel, name);
    } else if (name == "NAD-adapt") {
        return NAD_Adapt_Kernel(kmodel, "NAD");
    } else if (name == "NAD-adaptM") {
        return NAD_AdaptM_Kernel(kmodel, "NAD");
    } else {
        throw std::runtime_error("unknown solver name");
    }
    return nullptr;
}

std::shared_ptr<Solver> defaultSolverFactory(const std::string& name, std::shared_ptr<System> sys) {
    std::shared_ptr<Solver> ksolver = defaultSolverFactory(name, sys->getModel());
    ksolver->getSolverKernel()->setInputParam(sys->getParam());
    ksolver->getSolverKernel()->setInputDataSet(sys->getDataSet());
    return ksolver;
}

};  // namespace PROJECT_NS

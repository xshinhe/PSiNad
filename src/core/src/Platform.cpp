#include "kids/Platform.h"

#include "version.h"

namespace PROJECT_NS {

Platform::Platform(const std::string& name) : _platname{name} { _version = repo_version; }

std::shared_ptr<Platform> Platform::usePlatform(const std::string& platform_name) {
    if (false) {
    } else if (platform_name == "CPU") {
        return std::shared_ptr<CPU_Platform>(new CPU_Platform());
    } else if (platform_name == "CPU-MPI") {
        return std::shared_ptr<MPI_Platform>(new MPI_Platform());
    } else if (platform_name == "GPU") {
        return std::shared_ptr<GPU_Platform>(new GPU_Platform());
    } else if (platform_name == "CUDA") {
        return std::shared_ptr<CUDA_Platform>(new CUDA_Platform());
    } else if (platform_name == "OPENCL") {
        return std::shared_ptr<OPENCL_Platform>(new OPENCL_Platform());
    } else if (platform_name == "FPGA") {
        return std::shared_ptr<FPGA_Platform>(new FPGA_Platform());
    } else if (platform_name == "SW9") {
        return std::shared_ptr<SW9_Platform>(new SW9_Platform());
    }
    return std::shared_ptr<Platform>(new Platform("Unknown"));
}
};  // namespace PROJECT_NS

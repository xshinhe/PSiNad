#ifndef KIDS_PLATFORM_H
#define KIDS_PLATFORM_H

#include <memory>
#include <string>

namespace PROJECT_NS {

class Platform {
   public:
    inline std::string getName() { return _platname; }

    static std::shared_ptr<Platform> usePlatform(const std::string& platform_name);

   protected:
    Platform(const std::string& name);

    // virtual void initialize();

    // virtual void finalize();

    std::string _platname;

    std::string _version;
};

class CPU_Platform final : public Platform {
   public:
    CPU_Platform() : Platform("CPU"){};
};

class MPI_Platform final : public Platform {
   public:
    MPI_Platform() : Platform("CPU-MPI"){};
};

class GPU_Platform final : public Platform {
   public:
    GPU_Platform() : Platform("GPU"){};
};

class CUDA_Platform final : public Platform {
   public:
    CUDA_Platform() : Platform("CUDA"){};
};

class OPENCL_Platform final : public Platform {
   public:
    OPENCL_Platform() : Platform("OPENCL"){};
};

class FPGA_Platform final : public Platform {
   public:
    FPGA_Platform() : Platform("FPGA"){};
};

class SW9_Platform final : public Platform {
   public:
    SW9_Platform() : Platform("SW9"){};
};

// for python, cannot use c++'s MPI directly. One should inherit this class and override the functions

};  // namespace PROJECT_NS

#endif  // KIDS_PLATFORM_H

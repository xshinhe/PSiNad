#include "kids/Kernel.h"
#include "kids/Kernel_Elec_Functions.h"
#include "kids/Kernel_Load_DataSet.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Read_Dimensions.h"
#include "kids/Model.h"
#include "kids/Sampling_Elec.h"
#include "kids/Sampling_Nucl.h"
#include "kids/Solver.h"

namespace PROJECT_NS {

std::shared_ptr<Solver> Sampling_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name) {
    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(Kernel_name));

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("BAOAB_Integrator"));

    ker->appendChild(kmodel);

    std::shared_ptr<Solver> sol(new Solver(ker));
    return sol;
}

};  // namespace PROJECT_NS

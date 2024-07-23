#include "kids/Kernel.h"
#include "kids/Kernel_Dump_DataSet.h"
#include "kids/Kernel_Elec_Functions.h"
#include "kids/Kernel_Load_DataSet.h"
#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Read_Dimensions.h"
#include "kids/Kernel_Representation.h"
#include "kids/Model.h"
#include "kids/Sampling_Elec.h"
#include "kids/Sampling_Nucl.h"
#include "kids/Solver.h"

namespace PROJECT_NS {

std::shared_ptr<Solver> Sampling_Kernel(std::shared_ptr<Model> kmodel, std::string Kernel_name) {
    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(Kernel_name));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_NAForce>        knaf(new Kernel_NAForce());
    std::shared_ptr<Kernel_Elec_Functions> kfuncs(new Kernel_Elec_Functions());

    ker->appendChild(std::shared_ptr<Kernel_Load_DataSet>(new Kernel_Load_DataSet()))
        .appendChild(std::shared_ptr<Kernel_Random>(new Kernel_Random()))
        .appendChild(std::shared_ptr<Kernel_Read_Dimensions>(new Kernel_Read_Dimensions()))
        .appendChild(std::shared_ptr<Sampling_Nucl>(new Sampling_Nucl()))
        .appendChild(kmodel)
        .appendChild(krepr)
        .appendChild(kfuncs)
        .appendChild(knaf)
        .appendChild(std::shared_ptr<Sampling_Elec>(new Sampling_Elec()))
        .appendChild(std::shared_ptr<Kernel_Dump_DataSet>(new Kernel_Dump_DataSet()));
    std::shared_ptr<Solver> sol(new Solver(ker));
    return sol;
}

};  // namespace PROJECT_NS

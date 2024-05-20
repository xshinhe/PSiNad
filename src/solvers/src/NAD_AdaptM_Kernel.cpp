#include "kids/Kernel.h"
#include "kids/Kernel_Conserve.h"
#include "kids/Kernel_Dump_DataSet.h"
#include "kids/Kernel_Elec_NAD.h"
#include "kids/Kernel_Elec_Switch.h"
#include "kids/Kernel_GWP.h"
#include "kids/Kernel_Load_DataSet.h"
#include "kids/Kernel_NADForce.h"
#include "kids/Kernel_Prioritization.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Read_Dimensions.h"
#include "kids/Kernel_Recorder.h"
#include "kids/Kernel_Representation.h"
#include "kids/Kernel_Update.h"

namespace PROJECT_NS {

std::shared_ptr<Kernel> NAD_AdaptM_Kernel(std::shared_ptr<Kernel> kmodel, std::string NAD_Kernel_name) {
    bool take_ownership_false = false;

    int split = 4;

    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(NAD_Kernel_name));

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("MULTI_Integrator"));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_NADForce>       kforc(new Kernel_NADForce());

    std::shared_ptr<Kernel_Update_p> ku_p(new Kernel_Update_p(0.5e0 / (double) split));
    std::shared_ptr<Kernel_Update_x> ku_x(new Kernel_Update_x(0.5e0));
    std::shared_ptr<Kernel_Update_c> ku_c(new Kernel_Update_c(0.5e0 / (double) split));
    std::shared_ptr<Kernel>          kele;
    kele = std::shared_ptr<Kernel_Elec_NAD>(new Kernel_Elec_NAD(0.5e0 / (double) split));
    /// Result & Sampling & TCF
    std::shared_ptr<Kernel_Recorder> krecd(new Kernel_Recorder());

    for (int i = 0; i < split; ++i) {
        kinte->appendChild(ku_p);
        kinte->appendChild(krepr);
        kinte->appendChild(ku_c);
        kinte->appendChild(kele);
        kinte->appendChild(kforc);
    }

    kinte->appendChild(ku_x);
    kinte->appendChild(ku_x);
    kinte->appendChild(kmodel);
    kinte->appendChild(krepr);

    for (int i = 0; i < split; ++i) {
        kinte->appendChild(ku_c);
        kinte->appendChild(kele);
        kinte->appendChild(kforc);
        kinte->appendChild(ku_p);
        kinte->appendChild(krepr);
    }
    kinte->appendChild(std::shared_ptr<Kernel_Elec_Switch>(new Kernel_Elec_Switch()));
    kinte->appendChild(kforc);

    kinte->appendChild(std::shared_ptr<Kernel_Conserve>(new Kernel_Conserve()));

    std::shared_ptr<Kernel_Iter_Adapt> kiter(new Kernel_Iter_Adapt());
    kiter->appendChild(krecd);  // stacked in iteration
    kiter->appendChild(kinte);  // stacked in iteration

    // /// CMM kernel
    ker->appendChild(std::shared_ptr<Kernel_Load_DataSet>(new Kernel_Load_DataSet()))
        .appendChild(std::shared_ptr<Kernel_Random>(new Kernel_Random()))
        .appendChild(std::shared_ptr<Kernel_Read_Dimensions>(new Kernel_Read_Dimensions()))
        .appendChild(std::shared_ptr<Kernel_Prioritization>(new Kernel_Prioritization({kmodel, kinte}, 1)))
        .appendChild(std::shared_ptr<Kernel_Prioritization>(  //
            new Kernel_Prioritization({kmodel, krepr, kele, kforc, krecd}, 2)))
        .appendChild(kiter)
        .appendChild(std::shared_ptr<Kernel_Dump_DataSet>(new Kernel_Dump_DataSet()));
    return ker;
}

};  // namespace PROJECT_NS

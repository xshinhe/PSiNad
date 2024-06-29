#include "kids/Kernel.h"
#include "kids/Kernel_Conserve.h"
#include "kids/Kernel_Dump_DataSet.h"
#include "kids/Kernel_Elec_Functions.h"
#include "kids/Kernel_Elec_Switch.h"
#include "kids/Kernel_Load_DataSet.h"
#include "kids/Kernel_NAForce.h"
#include "kids/Kernel_Prioritization.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Read_Dimensions.h"
#include "kids/Kernel_Recorder.h"
#include "kids/Kernel_Representation.h"
#include "kids/Kernel_Update.h"
#include "kids/Model.h"
#include "kids/Solver.h"

namespace PROJECT_NS {

// CMM Solver Builder
std::shared_ptr<Solver> NAD_Kernel(std::shared_ptr<Model> kmodel, std::string NAD_Kernel_name) {
    bool take_ownership_false = false;

    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(NAD_Kernel_name));

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("BAOAB_Integrator"));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_Elec_Switch>    kswitch(new Kernel_Elec_Switch());
    std::shared_ptr<Kernel_NAForce>        knaf(new Kernel_NAForce());
    std::shared_ptr<Kernel_Elec_Functions> kfuncs(new Kernel_Elec_Functions());

    std::shared_ptr<Kernel_Update_p> ku_p(new Kernel_Update_p(0.5));
    std::shared_ptr<Kernel_Update_x> ku_x(new Kernel_Update_x(0.5));
    std::shared_ptr<Kernel_Update_U> ku_U(new Kernel_Update_U(1.0));

    /// Result & Sampling & TCF
    std::shared_ptr<Kernel_Recorder> krecd(new Kernel_Recorder());

    kinte->appendChild(ku_p);
    kinte->appendChild(ku_x);
    kinte->appendChild(ku_x);
    kinte->appendChild(kmodel);
    kinte->appendChild(krepr);
    kinte->appendChild(ku_U);
    kinte->appendChild(kswitch);

    /// kswarm = std::shared_ptr<Kernel_GWP>(new Kernel_GWP(kmodel, krepr, knaf));
    /// kinte->appendChild(kswarm);

    kinte->appendChild(knaf);
    kinte->appendChild(ku_p);
    kinte->appendChild(std::shared_ptr<Kernel_Conserve>(new Kernel_Conserve()));

    std::shared_ptr<Kernel_Iterative>   kiter(new Kernel_Iterative());
    std::shared_ptr<Kernel_Conditional> kcond(new Kernel_Conditional());
    kcond->appendChild(kfuncs);
    kcond->appendChild(krecd);
    kiter->appendChild(kcond);  // stacked in iteration
    kiter->appendChild(kinte);  // stacked in iteration

    // /// CMM kernel
    ker->appendChild(std::shared_ptr<Kernel_Load_DataSet>(new Kernel_Load_DataSet()))
        .appendChild(std::shared_ptr<Kernel_Random>(new Kernel_Random()))
        .appendChild(std::shared_ptr<Kernel_Read_Dimensions>(new Kernel_Read_Dimensions()))
        .appendChild(std::shared_ptr<Kernel_Prioritization>(new Kernel_Prioritization({kmodel, kinte}, 1)))
        .appendChild(std::shared_ptr<Kernel_Prioritization>(  //
            new Kernel_Prioritization({kmodel, krepr, kfuncs, knaf, krecd}, 2)))
        .appendChild(kiter)
        .appendChild(std::shared_ptr<Kernel_Dump_DataSet>(new Kernel_Dump_DataSet()));

    std::shared_ptr<Solver> sol(new Solver(ker));
    return sol;
}

};  // namespace PROJECT_NS

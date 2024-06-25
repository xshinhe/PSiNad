#include "kids/Kernel.h"
#include "kids/Kernel_Conserve.h"
#include "kids/Kernel_Dump_DataSet.h"
#include "kids/Kernel_Elec_NAD.h"
#include "kids/Kernel_GWP.h"
#include "kids/Kernel_Load_DataSet.h"
#include "kids/Kernel_NADForce.h"
#include "kids/Kernel_Prioritization.h"
#include "kids/Kernel_Random.h"
#include "kids/Kernel_Read_Dimensions.h"
#include "kids/Kernel_Recorder.h"
#include "kids/Kernel_Representation.h"
#include "kids/Kernel_Update.h"
#include "kids/Model.h"
#include "kids/Solver.h"

namespace PROJECT_NS {

std::shared_ptr<Solver> NAD_Adapt_Kernel(std::shared_ptr<Model> kmodel, std::string NAD_Kernel_name) {
    bool take_ownership_false = false;

    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(NAD_Kernel_name));

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("BAOAB_Integrator"));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_NADForce>       kforc(new Kernel_NADForce());

    std::shared_ptr<Kernel_Update_p> ku_p(new Kernel_Update_p(0.5));
    std::shared_ptr<Kernel_Update_x> ku_x(new Kernel_Update_x(0.5));
    std::shared_ptr<Kernel_Update_c> ku_c(new Kernel_Update_c(1.0));

    /// Result & Sampling & TCF
    std::shared_ptr<Kernel_Recorder> krecd(new Kernel_Recorder());

    kinte->appendChild(ku_p);
    kinte->appendChild(ku_x);
    kinte->appendChild(ku_x);
    kinte->appendChild(kmodel);
    kinte->appendChild(krepr);
    kinte->appendChild(ku_c);

    std::shared_ptr<Kernel> kele;
    if (false) {
        // } else if (NAD_Kernel_name == "CMM") {
        //     kele = std::shared_ptr<Kernel_Elec_CMM>(new Kernel_Elec_CMM());
        // } else if (NAD_Kernel_name == "SQC") {
        //     kele = std::shared_ptr<Kernel_Elec_SQC>(new Kernel_Elec_SQC());
        // } else if (NAD_Kernel_name == "MMD") {
        //     kele = std::shared_ptr<Kernel_Elec_MMD>(new Kernel_Elec_MMD());
        // } else if (NAD_Kernel_name == "SH") {
        //     kele = std::shared_ptr<Kernel_Hopping>(new Kernel_Hopping());
        // } else if (NAD_Kernel_name == "MMSH") {
        //     kele = std::shared_ptr<Kernel_Elec_MMSH>(new Kernel_Elec_MMSH());
    } else if (NAD_Kernel_name == "NAD") {
        kele = std::shared_ptr<Kernel_Elec_NAD>(new Kernel_Elec_NAD());
        // } else if (NAD_Kernel_name == "MCE") {
        //     kele = std::shared_ptr<Kernel_GWP>(new Kernel_GWP(kmodel, krepr, kforc));
    } else {
        throw std::runtime_error("unknown Elec Kernel");
    }
    kinte->appendChild(kele);

    kinte->appendChild(kforc);
    kinte->appendChild(ku_p);
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
    // return ker;

    std::shared_ptr<Solver> sol(new Solver(ker));
    return sol;
}

};  // namespace PROJECT_NS

#include "../kernels/Kernel_Conserve.h"
#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Dump_DataSet.h"
#include "../kernels/Kernel_Elec_CMM.h"
#include "../kernels/Kernel_Elec_CMSH.h"
#include "../kernels/Kernel_Elec_MMD.h"
#include "../kernels/Kernel_Elec_MMSH.h"
#include "../kernels/Kernel_Elec_SH.h"
#include "../kernels/Kernel_Elec_SQC.h"
#include "../kernels/Kernel_GWP.h"
#include "../kernels/Kernel_Load_DataSet.h"
#include "../kernels/Kernel_NADForce.h"
#include "../kernels/Kernel_Random.h"
#include "../kernels/Kernel_Record.h"
#include "../kernels/Kernel_Representation.h"
#include "../kernels/Kernel_Update.h"

namespace PROJECT_NS {

// CMM Solver Builder
std::shared_ptr<Kernel> NAD_Kernel(std::shared_ptr<Kernel> kmodel, std::string NAD_Kernel_name) {
    bool take_ownership_false = false;

    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel(NAD_Kernel_name));

    // Timer
    std::shared_ptr<Kernel_Timer> ktime(new Kernel_Timer());

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("BAOAB_Integrator"));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_NADForce> kforc(new Kernel_NADForce());

    std::shared_ptr<Kernel_Update_p> ku_p(new Kernel_Update_p(0.5));
    std::shared_ptr<Kernel_Update_x> ku_x(new Kernel_Update_x(0.5));
    std::shared_ptr<Kernel_Update_c> ku_c(new Kernel_Update_c(1.0));

    /// Result & Sampling & TCF
    std::shared_ptr<Kernel_Record> krecd(new Kernel_Record());

    kinte->push(ku_p);
    kinte->push(ku_x);
    kinte->push(ku_x);
    kinte->push(kmodel);
    kinte->push(krepr);
    kinte->push(ku_c);


    std::shared_ptr<Kernel> kele;
    if (false) {
    } else if (NAD_Kernel_name == "CMM") {
        kele = std::shared_ptr<Kernel_Elec_CMM>(new Kernel_Elec_CMM());
    } else if (NAD_Kernel_name == "SQC") {
        kele = std::shared_ptr<Kernel_Elec_SQC>(new Kernel_Elec_SQC());
    } else if (NAD_Kernel_name == "MMD") {
        kele = std::shared_ptr<Kernel_Elec_MMD>(new Kernel_Elec_MMD());
    } else if (NAD_Kernel_name == "SH") {
        kele = std::shared_ptr<Kernel_Elec_SH>(new Kernel_Elec_SH());
    } else if (NAD_Kernel_name == "MMSH") {
        kele = std::shared_ptr<Kernel_Elec_MMSH>(new Kernel_Elec_MMSH());
    } else if (NAD_Kernel_name == "CMSH") {
        kele = std::shared_ptr<Kernel_Elec_CMSH>(new Kernel_Elec_CMSH());
    } else if (NAD_Kernel_name == "MCE") {
        kele = std::shared_ptr<Kernel_GWP>(new Kernel_GWP(kmodel, krepr, kforc));
    } else {
        throw std::runtime_error("unknown Elec Kernel");
    }
    kinte->push(kele);

    kinte->push(kforc);
    kinte->push(ku_p);
    kinte->push(std::shared_ptr<Kernel_Conserve>(new Kernel_Conserve()));

    kinte->push(ktime);

    std::shared_ptr<Kernel_Iter> kiter(new Kernel_Iter());
    kiter->push(krecd);  // stacked in iteration
    kiter->push(kinte);  // stacked in iteration

    // /// CMM kernel
    ker->push(std::shared_ptr<Kernel_Load_DataSet>(new Kernel_Load_DataSet()))
        .push(std::shared_ptr<Kernel_Random>(new Kernel_Random()))
        .push(std::shared_ptr<Kernel_Declare>(new Kernel_Declare({kmodel, ktime, kinte})))
        .push(std::shared_ptr<Kernel_Initialize>(new Kernel_Initialize({kmodel, krepr, kele, krecd})))
        .push(kiter)
        .push(std::shared_ptr<Kernel_Dump_DataSet>(new Kernel_Dump_DataSet()));
    return ker;
}

};  // namespace PROJECT_NS

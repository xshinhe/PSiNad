#include "../kernels/Kernel_DataSetHandles.h"
#include "../kernels/Kernel_Dimension.h"
#include "../kernels/Kernel_Elec_SQC.h"
#include "../kernels/Kernel_NADForce.h"
#include "../kernels/Kernel_Random.h"
#include "../kernels/Kernel_Representation.h"
#include "../kernels/Kernel_Update.h"
#include "../kernels/some_Kernels.h"

namespace PROJECT_NS {

// CMM Solver Builder
std::shared_ptr<Kernel> SQC_Kernel(std::shared_ptr<Kernel> kmodel) {
    bool take_ownership_false = false;

    // Root Kernel
    std::shared_ptr<Kernel> ker(new Kernel("SQC"));
    std::shared_ptr<Kernel_Elec_SQC> kele(new Kernel_Elec_SQC());

    // Timer
    std::shared_ptr<Kernel_Timer> ktime(new Kernel_Timer());

    /// Integrator Kernel
    std::shared_ptr<Kernel> kinte(new Kernel("BAOAB_Integrator"));

    std::shared_ptr<Kernel_Representation> krepr(new Kernel_Representation());
    std::shared_ptr<Kernel_NADForce> kforc(new Kernel_NADForce());

    std::shared_ptr<Kernel_Update_p> ku_p(new Kernel_Update_p(0.5));
    std::shared_ptr<Kernel_Update_x> ku_x(new Kernel_Update_x(0.5));
    std::shared_ptr<Kernel_Update_rho> ku_rho(new Kernel_Update_rho(1.0));

    /// Result & Sampling & TCF
    std::shared_ptr<Kernel_Record> krec(new Kernel_Record());

    kinte->push(ku_p);
    kinte->push(ku_x);
    kinte->push(ku_x);
    kinte->push(kmodel);
    kinte->push(krepr);
    kinte->push(ku_rho);
    kinte->push(kele);
    kinte->push(kforc);
    kinte->push(ku_p);
    kinte->push(ktime);

    std::shared_ptr<Kernel_Iter> kiter(new Kernel_Iter());
    kiter->push(krec);   // stacked in iteration
    kiter->push(kinte);  // stacked in iteration

    // /// CMM kernel
    ker->push(std::shared_ptr<Kernel_Load_DataSet>(new Kernel_Load_DataSet()))
        .push(std::shared_ptr<Kernel_Random>(new Kernel_Random()))
        .push(std::shared_ptr<Kernel_Dimension>(new Kernel_Dimension()))
        .push(std::shared_ptr<Kernel_Declare>(new Kernel_Declare({kmodel, ktime, kinte})))
        .push(kiter)
        .push(std::shared_ptr<Kernel_Dump_DataSet>(new Kernel_Dump_DataSet()));
    return ker;
}

};  // namespace PROJECT_NS

#include "psnd/Kernel.h"
#include "psnd/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(SwarmCoupPolicy,
              hcps,  // hydrodynamics cps
              ccps   // coupled-cps
);                   // Read from dataset

class Kernel_Update_Coup final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_Update_Coup(double scale) : scale{scale} {};

   private:
    psnd_real sigma_nuc, sigma_ele;

    SwarmCoupPolicy::_type swarm_type;

    psnd_real scale;
    psnd_real dt;

    span<psnd_real> relwgt, relwgt2;
    span<psnd_real> gf_x, gf_p, gf_c, gf_all, avgx, varx, avgp, varp, avgxf, varxf;
    span<psnd_real> xintercept, xinterceptf, xslope;
    span<psnd_real> term_1, term_2, fadiat, pb;

    span<psnd_real>    x, p, m, f;
    span<psnd_real>    dV, dE, T, T_init;
    span<psnd_complex> c, rho_ele, U, Ucdt;
    span<psnd_complex> K1, K2, rhored;


    double xi1, xi2, gamma1, gamma2;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

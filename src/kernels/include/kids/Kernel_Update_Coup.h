#include "kids/Kernel.h"
#include "kids/Policy.h"

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
    kids_real sigma_nuc, sigma_ele;

    SwarmCoupPolicy::_type swarm_type;

    kids_real scale;
    kids_real dt;

    span<kids_real> relwgt, gf_x, gf_p, gf_c, avgx, varx, avgp, varp, avgxf, varxf;
    span<kids_real> xintercept, xinterceptf, xslope;
    span<kids_real> term_1, term_2, fadiat, pb;

    span<kids_real>    x, p, m, f;
    span<kids_real>    dV, dE, T, T_init;
    span<kids_complex> c, rho_ele, U, Udt;
    span<kids_complex> K1, K2, rhored;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

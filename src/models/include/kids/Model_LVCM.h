#ifndef Model_LVCM_H
#define Model_LVCM_H

#include "kids/Model.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(LVCMPolicy,  //
              PYR3,        // 3-mode pyrazine model
              PYR3_SPEC,   // 3-mode pyrazine model spectrum
              PYR4,        // 4-mode pyrazine model
              PYR4_SPEC,   // 4-mode pyrazine model spectrum
              PYR24,       // 24-mode pyrazine model
              CRC2,        // 2-mode Cr(CO)5 model
              CRC5,        // 5-mode Cr(CO)5 model
              BEN5,        // 5-mode Benzene model
              BUTA5,       // 5-mode butatiene model ??
              PENTA5,      // 5-mode pentaxxxxene model
              CED2,        // 2-state atom-in-cavity model
              CED3,        // 3-state atom-in-cavity model
              PYR2CED,     // pyrazine-in-cavity model
              Read);       //

class Model_LVCM final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_LVCM(){};

   private:
    LVCMPolicy::_type lvcm_type;
    int               N_coup;
    int               N_mode;

    span<kids_real> Hsys;
    span<kids_real> kcoeff, lcoeff;

    span<kids_real> x_sigma;
    span<kids_real> p_sigma;
    span<kids_real> x0;
    span<kids_real> p0;

    // integrator
    span<kids_real> x, p, m, w;

    // model
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess;
    span<kids_real> V, dV, ddV;

    span<kids_real> Kmat, Qmat, Tmod;

    // int N_ligh;
    // N = N_mode + N_coup + N_ligh

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_LVCM_H

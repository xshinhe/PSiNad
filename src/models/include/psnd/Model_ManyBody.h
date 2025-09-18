#ifndef Model_ManyBody_H
#define Model_ManyBody_H

#include "psnd/Model.h"
#include "psnd/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(ManyBodyPolicy,  //
              TwoSite,         // two-body composite mode
              Heisenberg,      // Heisenberg model
              Ising,           // Ising model
              XY,              // XY model
              Hubbard,         // Hubbard model
              Anderson,        // Anderson model
              BowTie,          // Bow-Tie model
              DemkovOsherov,   // Demkov-Osherov model
              TavisCummings,   // Tavis-Cummings model -> Jaynes-Cummings model
              Dicke,           // Dicke model -> Rabi model
              Read);           // (read from file)

class Model_ManyBody : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_ManyBody() {};

   private:
    span<psnd_real> w;
    span<psnd_real> x_sigma;
    span<psnd_real> p_sigma;
    span<psnd_real> x0;
    span<psnd_real> p0;
    span<psnd_real> x, p;
    span<psnd_real> mass;
    span<psnd_real> vpes, grad, hess;
    span<psnd_real> V, dV, ddV;

    span<psnd_real>    Jpmat;
    span<psnd_real>    Jzmat;
    span<psnd_complex> SXred;
    span<psnd_complex> SYred;
    span<psnd_complex> SZred;
    span<psnd_complex> H1, H2;
    span<psnd_complex> H;
    span<psnd_complex> rho_ele;

    psnd_real Jp, Jz;
    psnd_real alpha;
    psnd_real omega;

    ManyBodyPolicy::_type manybody_type;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& executeKernel_impl(Status& stat);

    Status& execute_Heisenberg(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_ManyBody_H
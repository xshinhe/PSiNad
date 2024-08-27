#ifndef Model_ManyBody_H
#define Model_ManyBody_H

#include "kids/Model.h"
#include "kids/Policy.h"

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

    Model_ManyBody(){};

   private:
    span<kids_real> w;
    span<kids_real> x_sigma;
    span<kids_real> p_sigma;
    span<kids_real> x0;
    span<kids_real> p0;
    span<kids_real> x, p;
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess;
    span<kids_real> V, dV, ddV;

    span<kids_real>    Jpmat;
    span<kids_real>    Jzmat;
    span<kids_complex> SXred;
    span<kids_complex> SYred;
    span<kids_complex> SZred;
    span<kids_complex> H1, H2;
    span<kids_complex> H;
    span<kids_complex> rho_ele;

    kids_real Jp, Jz;
    kids_real alpha;
    kids_real omega;

    ManyBodyPolicy::_type manybody_type;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& executeKernel_impl(Status& stat);

    Status& execute_Heisenberg(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Model_ManyBody_H
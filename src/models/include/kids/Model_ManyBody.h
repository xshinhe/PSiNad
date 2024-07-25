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
    kids_real* w;
    kids_real* x_sigma;
    kids_real* p_sigma;
    kids_real* x0;
    kids_real* p0;
    kids_real *x, *p;
    kids_real* mass;
    kids_real *vpes, *grad, *hess;
    kids_real *V, *dV, *ddV;

    kids_real*    Jpmat;
    kids_real*    Jzmat;
    kids_complex* SXred;
    kids_complex* SYred;
    kids_complex* SZred;
    kids_complex *H1, *H2;
    kids_complex* H;
    kids_complex* rho_ele;

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
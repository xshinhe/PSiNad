#ifndef Model_HarmonicBath_H
#define Model_HarmonicBath_H

#include "kids/Model.h"
#include "kids/Policy.h"


namespace PROJECT_NS {


DEFINE_POLICY(HarmonicBathPolicy,  //
              Debye,               //
              Ohmic,               //
              OhmicET,             //
              Closure,             //
              HuangRhys,           //
              GFactor,             //
              ReadFormula,         //
              ReadFile,            //
              Read);               //

DEFINE_POLICY(StrengthPolicy,  //
              Lambda,          //
              Alpha,           //
              Eta,             //
              Erg);            //

class Model_HarmonicBath final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    bool is_classical = false;  // for classical / quantum fluctuation

    bool is_correlated = false;  // for normal-mode / not-normal-mode treatment

    bool is_et_transform = false;  // for ET bath model

    int                       nbath;
    int                       Nb;
    HarmonicBathPolicy::_type bath_type;
    StrengthPolicy::_type     strength_type;
    double                    omegac;
    double                    lambda;
    double                    beta;

    double J_Debye(double w);

    double J_Ohmic(double w);

    static double J(double w, double* w_arr = nullptr, double* c_arr = nullptr, int Nb = 0);

    static int fun_Cw(kids_complex* Cw_arr, double* w, int Nw, double* w_arr, double* c_arr, double beta, int Nb);

   private:
    kids_real* Kmat;
    kids_real* Tmod;
    kids_real* w;

    kids_real* coeffs;
    kids_real* omegas;
    kids_real *x_sigma, *x0;
    kids_real *p_sigma, *p0;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);
};


};  // namespace PROJECT_NS

#endif  // Model_HarmonicBath_H
#ifndef LVCM_H
#define LVCM_H

#include "../forcefieldbase.h"

namespace PROJECT_NS {


namespace LVCMPolicy {
enum _enum { PYR3, PYR3_SPECTRUM, PYR4, PYR4_SPECTRUM, PYR24, PYR24_SPECTRUM, CRC2, CRC5, BEN5 };
const std::map<std::string, _enum> _dict = {
    {"pyr3", PYR3},
    {"pyr3_spectrum", PYR3_SPECTRUM},
    {"pyr4", PYR4},
    {"pyr4_spectrum", PYR4_SPECTRUM},
    {"pyr24", PYR24_SPECTRUM},
    {"crc2", CRC2},
    {"crc5", CRC5},
    {"ben5", BEN5},
};
};  // namespace LVCMPolicy

class LVCM_ForceField final : public Model {
   public:
    static inline std::string name() { return "lvcm"; }

   private:
    int N, F, NN, FF, NNF, NFF, NNFF;

    // external variables
    double *x, *p, *f;

    // shared variables
    double *m, *w, *x_zero, *p_zero;
    double *vpes, *grad, *hess;
    double *V, *dV, *ddV;

    // internal variables
    double* Hsys;


    Option ff_option;

    virtual void read_param_impl(Param* P) {
        N    = P->get<int>("N", LOC());   // get DOFs
        F    = P->get<int>("F", LOC());   // get DOFs
        Nb   = P->get<int>("Nb", LOC());  // get DOFs
        NN   = N * N;
        FF   = F * F;
        NNF  = NN * F;
        NFF  = N * FF;
        NNFF = NN * FF;

        ff_option.flag = P->get<std::string>("ff", LOC(), "ced2");
        ff_option.type = AtomCED_ff_policy::_dict.at(ff_option.flag);
    }

    virtual void init_data_impl(State* S) {
        // external variables
        x = S->reg<double>("integrator.x", N);
        p = S->reg<double>("integrator.p", N);
        f = S->reg<double>("integrator.f", N);

        // shared variables
        vpes = S->reg<double>("model.vpes");
        grad = S->reg<double>("model.grad", N);
        hess = S->reg<double>("model.hess", NN);
        V    = S->reg<double>("model.V", FF);
        dV   = S->reg<double>("model.dV", NFF);
        // ddV    = S->reg<double>("model.ddV", NNFF);
        m      = S->reg<double>("model.m", N);
        w      = S->reg<double>("model.w", N);
        x_zero = S->reg<double>("model.x_zero", N);
        p_zero = S->reg<double>("model.p_zero", N);

        // internal variables
        Hsys   = S->reg<double>("model.hidden.Hsys", FF);
        omegas = S->reg<double>("model.hidden.omegas", Nb);
        coeffs = S->reg<double>("model.hidden.coeffs", Nb);

        // initialization of hidden variables
        switch (ff_option.type) {
            case AtomCED_ff_policy::CED2:
                // CHECK_EQ(F, 2);
                Hsys[0] = -0.6738f;  // e1
                Hsys[1] = +1.034f;   // mu12
                Hsys[2] = +1.034f;   // mu21
                Hsys[3] = -0.2798f;  // e2
                Lcav    = 2.362e5;   // a.u.
                Rcav    = Lcav / 2;
                break;
            case AtomCED_ff_policy::CED3:
                // CHECK_EQ(F, 3);
                Hsys[0] = -0.6738f;  // e1
                Hsys[1] = +1.034f;   // mu12
                Hsys[2] = +0.000f;   // mu13
                Hsys[3] = +1.034f;   // mu21
                Hsys[4] = -0.2798f;  // e2
                Hsys[5] = -2.536f;   // mu23
                Hsys[6] = +0.000f;   // mu31
                Hsys[7] = -2.536f;   // mu32
                Hsys[8] = -0.1547f;  // e3
                Lcav    = 2.362e5;   // a.u.
                Rcav    = Lcav / 2;
                break;
            default: {
                // LOG(FATAL);
            }
        }
        for (int j = 0; j < Nb; ++j) {
            omegas[j] = (2 * j + 1) * lightspeed * phys::math::pi / Lcav;
            coeffs[j] = sqrt(2.0f / (epsilon0 * Lcav)) * sin((2 * j + 1) * phys::math::pi * Rcav / Lcav);
        }

        // initialization of shared variables
        for (int j = 0; j < N; ++j) {
            m[j]      = 1.0f;
            w[j]      = omegas[j];
            x_zero[j] = 0.0f;
            p_zero[j] = 0.0f;
        }
    }


   public:
    LVCM_ForceField(const Param& iparm, const int& child);
    LVCM_ForceField(const Param& iparm);
    LVCM_ForceField(const std::string& iparm_str) : LVCM_ForceField(Param::parse(iparm_str)){};

    virtual ~LVCM_ForceField(){};

    static inline std::string name() { return "lvcm"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);


    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);


    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    // virtual int ForceField_epes_PYR3(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const
    // int& rdim,
    //                                  const int& fdim);

    // virtual int ForceField_epes_PYR3_SPECTRUM(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
    //                                           const int& rdim, const int& fdim);

    virtual int ForceField_epes_PYR(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    // virtual int ForceField_epes_PYR_SPECTRUM(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
    //                                          const int& rdim, const int& fdim);

    virtual int ForceField_epes_CrCO5_2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_CrCO5_5(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_BEN_5(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                      const int& rdim, const int& fdim);


   protected:
    DEFINE_POINTER_PROTECTED(num_real, Hsys);
    DEFINE_POINTER_PROTECTED(num_real, ECI);
    DEFINE_POINTER_PROTECTED(num_real, KCI);
    double lcoeff;
    bool first_call = true;

    int fftype;
    std::string ffflag;
};

};  // namespace PROJECT_NS

#endif  // LVCM_H

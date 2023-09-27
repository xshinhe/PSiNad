#ifndef AtomCED_H
#define AtomCED_H

// #include "../forcefieldbase.h"
#include "../../core/Kernel.h"
#include "../../core/Param.h"
#include "../../core/State.h"
#include "../../kernels/some_Kernels.h"


namespace PROJECT_NS {

DEFINE_POLICY(AtomCED, CED2, CED3, PYR3)

namespace AtomCED_ff_policy {
enum _enum { CED2, CED3 };
const std::map<std::string, _enum> _dict = {{"ced2", CED2}, {"ced3", CED3}};
};  // namespace AtomCED_ff_policy

class AtomCED_ForceField final : public Kernel {
   public:
    static inline std::string name() { return "atomced"; }

   private:
    const double lightspeed = 137.03599907444f;
    const double epsilon0   = 0.25f / phys::math::pi;

    int N, F, NN, FF, NNF, NFF, NNFF;

    // external variables
    double *x, *p, *f;

    // shared variables
    double *m, *w, *x_zero, *p_zero;
    double *vpes, *grad, *hess;
    double *V, *dV, *ddV;

    // internal variables
    double* Hsys;
    double* omegas;
    double* coeffs;

    Option ff_option;

    double Lcav, Rcav;
    int Nb;
    bool first_call = true;

    virtual void read_param_impl(Param* P) {
        N    = P->get<int>("N", LOC());   // get DOFs
        F    = P->get<int>("F", LOC());   // get DOFs
        Nb   = P->get<int>("Nb", LOC());  // get DOFs
        NN   = N * N;
        FF   = F * F;
        NNF  = NN * F;
        NFF  = N * FF;
        NNFF = NN * FF;

        ff_option.flag = P->get<std::string>("ff", LOC(), "diabatic");
        ff_option.type = AtomCED_ff_policy::_dict.at(ff_option.flag);
    }

    virtual void init_data_impl(State* S) {
        // external variables
        x = S->reg<double>("integrator.x", N);
        p = S->reg<double>("integrator.p", N);
        f = S->reg<double>("integrator.f", N);

        // shared variables
        vpes   = S->reg<double>("model.vpes");
        grad   = S->reg<double>("model.grad", N);
        hess   = S->reg<double>("model.hess", NN);
        V      = S->reg<double>("model.V", FF);
        dV     = S->reg<double>("model.dV", NFF);
        ddV    = S->reg<double>("model.ddV", NNFF);
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

    virtual int exec_kernel_impl(int stat = -1) {
        // vpes
        vpes[0] = 0.0f;
        for (int j = 0; j < N; ++j) { vpes[0] += 0.5f * (p[j] * p[j] / m[j] + m[j] * m[j] * w[j] * x[j] * x[j]); }

        // grad
        for (int j = 0; j < N; ++j) { grad[j] = m[j] * w[j] * w[j] * x[j]; }

        // hess
        for (int j = 0; j < NN; ++j) hess[j] = 0;
        int add = N + 1;
        for (int j = 0, idx = 0; j < N; ++j, idx += add) hess[idx] = m[j] * w[j] * w[j];

        // V
        for (int i = 0, idxV = 0; i < F; ++i)
            for (int k = 0; k < F; ++k, ++idxV) V[idxV] = (i == k) ? Hsys[idxV] : 0.0f;
        for (int j = 0; j < N; ++j) {
            for (int i = 0, idxV = 0; i < F; ++i) {
                for (int k = 0; k < F; ++k, ++idxV) {
                    if (i == k) continue;
                    V[idxV] += omegas[j] * x[j] * coeffs[j] * Hsys[idxV];
                }
            }
        }

        // dV
        if (first_call) {
            // dV
            for (int j = 0, idxdV = 0; j < N; ++j) {
                for (int i = 0, idxV = 0; i < F; ++i) {
                    for (int k = 0; k < F; ++k, ++idxV, ++idxdV) {
                        dV[idxdV] = (i == k) ? 0.0f : omegas[j] * coeffs[j] * Hsys[idxV];
                    }
                }
            }
        }

        // ddV
        for (int i = 0; i < NNFF; ++i) ddV[i] = 0;

        first_call = false;
        return 0;
    };

    virtual Kernel* new_model_init_kernel() { return new Kernel_Init_Nucl_Wigner_Gaussian(); }
};

};  // namespace PROJECT_NS


#endif  // AtomCED_H

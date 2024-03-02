#ifndef MODEL_MD1D_H
#define MODEL_MD1D_H

#include "../core/Kernel.h"
#include "../kernels/Kernel_Random.h"

namespace kids {


class Model_HO final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_HO"; }

   private:
    int size;
    double *x, *p;
    double *x_init, *p_init;
    double *f;
    double *m, *w;

    double beta;

    virtual void read_param_impl(Param *PM) {
        size = PM->get<int>("N", LOC());             //
        beta = PM->get<double>("beta", LOC(), 1.0);  //
    }

    virtual void init_data_impl(DataSet *DS) {
        // external
        x = DS->def<double>("integrator.x", size);
        p = DS->def<double>("integrator.p", size);
        f = DS->def<double>("integrator.f", size);

        x_init = DS->def<double>("init.x", size);
        p_init = DS->def<double>("init.p", size);

        // shared
        w = DS->def<double>("model.w", size);
        m = DS->def<double>("model.m", size);
        for (int i = 0; i < size; ++i) m[i] = 1, w[i] = 1;
    }

    virtual void init_calc_impl(int stat = -1) {
        for (int i = 0; i < size; ++i) {
            double Qi_div_beta = 0.5f * w[i] / std::tanh(0.5f * beta * w[i]);
            double xi_sigma    = std::sqrt(Qi_div_beta / (m[i] * w[i] * w[i]));
            double pi_sigma    = std::sqrt(Qi_div_beta * m[i]);

            double u1, u2;
            Kernel_Random::rand_gaussian(&u1);
            Kernel_Random::rand_gaussian(&u2);
            x_init[i] = u1 * xi_sigma;
            p_init[i] = u2 * pi_sigma;
        }
        // copy to external
        for (int i = 0; i < size; ++i) {
            x[i] = x_init[i];
            p[i] = p_init[i];
        }
    }

    virtual int exec_kernel_impl(int stat = -1) {
        for (int i = 0; i < size; ++i) f[i] = m[i] * w[i] * w[i] * x[i];
        return 0;
    }
};

class MODEL_MD1D final : public Kernel {
   public:
    inline virtual const std::string name() { return "MODEL_MD1D"; }

   private:
    int size;
    double *x, *f;
    double *m, *w;

    virtual void read_param_impl(Param *PM) {
        size = PM->get<int>("N", LOC());  //
    }

    virtual void init_data_impl(DataSet *DS) {
        // external
        x = DS->def<double>("integrator.x", size);
        f = DS->def<double>("integrator.f", size);
        // shared
        w = DS->def<double>("model.w", size);
        m = DS->def<double>("model.m", size);
        for (int i = 0; i < size; ++i) m[i] = 1, w[i] = 1;
    }

    virtual int exec_kernel_impl(int stat = -1) {
        for (int i = 0; i < size; ++i) f[i] = m[i] * w[i] * w[i] * x[i];
        return 0;
    }
};
};  // namespace kids

#endif  // MODEL_MD1D_H

#ifndef Kernel_Update_H
#define Kernel_Update_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Declare final : public Kernel {
   public:
    Kernel_Declare(std::shared_ptr<Kernel> ker) : Kernel() { _ref_kernels.push_back(ker); }

    Kernel_Declare(std::vector<std::shared_ptr<Kernel>> kers) : Kernel() {
        for (auto& ker : kers) { _ref_kernels.push_back(ker); }
    }

    inline virtual const std::string name() {
        std::stringstream ss;
        ss << "Kernel_Declare";
        for (auto& ker : _ref_kernels) ss << " #" << std::setfill('0') << std::setw(2) << ker->id();
        return ss.str();
    }

   private:
    // std::shared_ptr<Kernel> _ref_kernel;
    std::vector<std::shared_ptr<Kernel>> _ref_kernels;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);
};

/**
 * @brief iterative kernel wrapper/(interface) for other kernels
 */
class Kernel_Iter final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Iter"; }

   private:
    int* istep_ptr;
    int* nstep_ptr;

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};

class Kernel_Timer final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Timer"; }

   private:
    double t, *t_ptr;
    double t0, tend, dt;
    int sstep, *sstep_ptr;
    int istep, *istep_ptr, nstep, *nstep_ptr;
    int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

class Kernel_Update_x final : public Kernel {
   public:
    Kernel_Update_x(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_x"; }

   private:
    double *x, *p, *m, *minv;

    double scale;
    double dt, sdt;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};

class Kernel_Update_p : public Kernel {
   public:
    Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_p"; }

   private:
    double *p, *f;

    double scale;
    double dt, sdt;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};


class Kernel_Update_T : public Kernel {
   public:
    Kernel_Update_T(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_T"; }

   private:
    double *p, *m;
    // for Langevin
    double *c1, *c2p;
    // for NHC
    double *nhc_x, *nhc_p, *nhc_G, *nhc_Q;
    //
    double scale;
    double dt, sdt;
    double beta;
    double gammal;
    double randu;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* S);

    virtual int exec_kernel_impl(int stat = -1);
};

class Kernel_Update_rho final : public Kernel {
   public:
    Kernel_Update_rho(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_rho"; }

   private:
    num_complex* c;
    num_complex* rho_ele;
    num_complex *U, *invexpiEdt, *invexpiLdt;

    num_real *V, *dV, *ddV;
    num_real *E, *T, *dE, *ddE;
    num_real* L;
    num_complex *R, *dL, *ddL;
    num_complex *H, *dH, *ddH;
    num_complex* Matr;

    double scale;
    double dt, sdt;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* S);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_H

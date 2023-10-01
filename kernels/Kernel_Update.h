#ifndef Kernel_Update_H
#define Kernel_Update_H

#include "Kernel_Iter.h"
#include "Kernel_Timer.h"
#include "Kernel_Update_T.h"
#include "Kernel_Update_c.h"
#include "Kernel_Update_p.h"
#include "Kernel_Update_x.h"

// #include "../core/Kernel.h"

// namespace PROJECT_NS {

// /**
//  * @brief iterative kernel wrapper/(interface) for other kernels
//  */
// class Kernel_Iter final : public Kernel {
//    public:
//     inline virtual const std::string name() { return "Kernel_Iter"; }

//    private:
//     int* istep_ptr;
//     int* nstep_ptr;

//     virtual void init_data_impl(DataSet* DS);

//     virtual int exec_kernel_impl(int stat = -1);
// };

// class Kernel_Timer final : public Kernel {
//    public:
//     inline virtual const std::string name() { return "Kernel_Timer"; }

//    private:
//     double t, *t_ptr;
//     double t0, tend, dt;
//     int sstep, *sstep_ptr;
//     int istep, *istep_ptr, nstep, *nstep_ptr;
//     int isamp, *isamp_ptr, nsamp, *nsamp_ptr;

//     virtual void read_param_impl(Param* PM);

//     virtual void init_data_impl(DataSet* DS);

//     virtual void init_calc_impl(int stat = -1);

//     virtual int exec_kernel_impl(int stat = -1);
// };

// class Kernel_Update_x final : public Kernel {
//    public:
//     Kernel_Update_x(double scale) : Kernel(), scale{scale} {};

//     inline virtual const std::string name() { return "Kernel_Update_x"; }

//    private:
//     double *x, *p, *m, *minv;

//     double scale;
//     double dt, sdt;

//     virtual void read_param_impl(Param* PM);

//     virtual void init_data_impl(DataSet* DS);

//     virtual int exec_kernel_impl(int stat = -1);
// };

// class Kernel_Update_p : public Kernel {
//    public:
//     Kernel_Update_p(double scale) : Kernel(), scale{scale} {};

//     inline virtual const std::string name() { return "Kernel_Update_p"; }

//    private:
//     double *p, *f;

//     double scale;
//     double dt, sdt;

//     virtual void read_param_impl(Param* PM);

//     virtual void init_data_impl(DataSet* DS);

//     virtual int exec_kernel_impl(int stat = -1);
// };


// class Kernel_Update_T : public Kernel {
//    public:
//     Kernel_Update_T(double scale) : Kernel(), scale{scale} {};

//     inline virtual const std::string name() { return "Kernel_Update_T"; }

//    private:
//     double *p, *m;
//     // for Langevin
//     double *c1, *c2p;
//     // for NHC
//     double *nhc_x, *nhc_p, *nhc_G, *nhc_Q;
//     //
//     double scale;
//     double dt, sdt;
//     double beta;
//     double gammal;
//     double randu;

//     virtual void read_param_impl(Param* PM);

//     virtual void init_data_impl(DataSet* S);

//     virtual int exec_kernel_impl(int stat = -1);
// };

// class Kernel_Update_c final : public Kernel {
//    public:
//     Kernel_Update_c(double scale) : Kernel(), scale{scale} {};

//     inline virtual const std::string name() { return "Kernel_Update_c"; }

//    private:
//     num_complex* Udt;  ///< short time propagator
//     num_complex* U;    ///< full propagator along classical path approximation (CPA)

//     ///< solve Diabatic propagator
//     num_real* E;  ///< Eigenvalue for diabatic V
//     num_real* T;  ///< Eigenvector for diabatic V

//     ///< solve Adiabatic propagator
//     num_real* L;     ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
//     num_complex* R;  ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

//     num_complex* invexpidiagdt;  ///< temporary variables

//     double scale;
//     double dt, sdt;

//     virtual void read_param_impl(Param* PM);

//     virtual void init_data_impl(DataSet* S);

//     virtual int exec_kernel_impl(int stat = -1);
// };


// };  // namespace PROJECT_NS


#endif  // Kernel_Update_H

#include "Kernel_Update.h"

#include "../core/linalg.h"
#include "Kernel_Dimension.h"
#include "Kernel_Elec.h"
#include "Kernel_Random.h"
#include "Kernel_Representation.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(4) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

namespace PROJECT_NS {

// std::shared_ptr<Kernel> _ref_kernel;
std::vector<std::shared_ptr<Kernel>> _ref_kernels;

void Kernel_Declare::read_param_impl(Param* PM) {
    for (auto& ker : _ref_kernels) ker->read_param(PM);
}

void Kernel_Declare::init_data_impl(DataSet* DS) {
    for (auto& ker : _ref_kernels) ker->init_data(DS);
}

void Kernel_Declare::init_calc_impl(int stat) {
    // for (auto& ker : _ref_kernels) ker->init_calc(stat);
}

void Kernel_Iter::init_data_impl(DataSet* DS) {
    istep_ptr = DS->reg<int>("timer.istep");
    nstep_ptr = DS->reg<int>("timer.nstep");
}

int Kernel_Iter::exec_kernel_impl(int stat) {
    int& istep_ref = *istep_ptr;
    int& nstep_ref = *nstep_ptr;
    while (istep_ref < nstep_ref) {
        for (auto& pkernel : _kernel_vector) { pkernel->exec_kernel(stat); }
    }
    return 0;
}

void Kernel_Timer::read_param_impl(Param* PM) {
    t0    = PM->get<double>("t0", LOC(), phys::time_d, 0.0f);
    tend  = PM->get<double>("tend", LOC(), phys::time_d, 1.0f);
    dt    = PM->get<double>("dt", LOC(), phys::time_d, 0.1f);
    sstep = PM->get<int>("sstep", LOC(), 1);
    nstep = sstep * (int((tend - t0) / (sstep * dt)) + 1);
    nsamp = nstep / sstep + 1;
}

void Kernel_Timer::init_data_impl(DataSet* DS) {
    t_ptr     = DS->reg<num_real>("integrator.t");
    sstep_ptr = DS->reg<int>("timer.sstep");
    istep_ptr = DS->reg<int>("timer.istep");
    nstep_ptr = DS->reg<int>("timer.nstep");
    isamp_ptr = DS->reg<int>("timer.isamp");
    nsamp_ptr = DS->reg<int>("timer.nsamp");
    // initial
    *sstep_ptr = sstep;
    *nstep_ptr = nstep;
    *nsamp_ptr = nsamp;
}

void Kernel_Timer::init_calc_impl(int stat) {
    (*istep_ptr) = 0;
    (*isamp_ptr) = 0;
    (*t_ptr)     = t0;
}

int Kernel_Timer::exec_kernel_impl(int stat) {
    (*t_ptr) += dt;
    (*istep_ptr)++;
    (*isamp_ptr) = (*istep_ptr) / sstep;
    return 0;
}

void Kernel_Update_x::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);
    sdt = scale * dt;
}

void Kernel_Update_x::init_data_impl(DataSet* DS) {
    x              = DS->reg<num_real>("integrator.x", Kernel_Dimension::N);
    p              = DS->reg<num_real>("integrator.p", Kernel_Dimension::N);
    m              = DS->reg<num_real>("integrator.m", Kernel_Dimension::N);
    minv           = DS->reg<num_real>("integrator.minv", Kernel_Dimension::N);
    num_real* mass = DS->reg<num_real>("model.mass", Kernel_Dimension::N);
    for (int i = 0; i < Kernel_Dimension::N; ++i) { m[i] = mass[i], minv[i] = 1 / m[i]; }
}

int Kernel_Update_x::exec_kernel_impl(int stat) {
    for (int i = 0; i < Kernel_Dimension::N; ++i) x[i] += p[i] * minv[i] * sdt;
    // ARRAY_SHOW(x, 1, Kernel_Dimension::N);
    return 0;
}



void Kernel_Update_p::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);  //
    sdt = scale * dt;
}

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    f = DS->reg<num_real>("integrator.f", Kernel_Dimension::N);
    p = DS->reg<num_real>("integrator.p", Kernel_Dimension::N);
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    for (int i = 0; i < Kernel_Dimension::N; ++i) p[i] -= f[i] * sdt;
    // ARRAY_SHOW(p, 1, Kernel_Dimension::N);
    return 0;
}



void Kernel_Update_T::read_param_impl(Param* PM) {
    dt     = PM->get<double>("dt", LOC(), phys::time_d);  //
    gammal = PM->get<double>("gammal", LOC(), 0.1);       //
    sdt    = scale * dt;
}

void Kernel_Update_T::init_data_impl(DataSet* S) {
    m = S->reg<num_real>("integrator.m", Kernel_Dimension::N);
    p = S->reg<num_real>("integrator.p", Kernel_Dimension::N);

    // if Langevin dynamics, set optimal c1 & c2p
    c1  = S->reg<num_real>("integrator.c1", Kernel_Dimension::N);
    c2p = S->reg<num_real>("integrator.c2p", Kernel_Dimension::N);
    for (int i = 0; i < Kernel_Dimension::N; ++i) {
        c1[i]  = exp(-gammal * dt);
        c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
    }

    // if for NHC; registeration for auxiliary variables
    // ...
}

int Kernel_Update_T::exec_kernel_impl(int stat) {
    for (int i = 0; i < Kernel_Dimension::N; ++i) {
        Kernel_Random::rand_gaussian(&randu);
        p[i] = c1[i] * p[i] + c2p[i] * sqrt(m[i] / beta) * randu;
    }
    return 0;
}


void Kernel_Update_rho::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);  //
    sdt = scale * dt;
}

void Kernel_Update_rho::init_data_impl(DataSet* S) {
    V          = S->reg<num_real>("model.V", Kernel_Dimension::FF);
    E          = S->reg<num_real>("model.rep.E", Kernel_Dimension::F);
    T          = S->reg<num_real>("model.rep.T", Kernel_Dimension::FF);
    L          = S->reg<num_real>("model.rep.L", Kernel_Dimension::F);
    R          = S->reg<num_complex>("model.rep.R", Kernel_Dimension::FF);
    H          = S->reg<num_complex>("model.rep.H", Kernel_Dimension::FF);
    U          = S->reg<num_complex>("integrator.U", Kernel_Dimension::FF);
    invexpiEdt = S->reg<num_complex>("integrator.invexpiEdt", Kernel_Dimension::F);
    invexpiLdt = S->reg<num_complex>("integrator.invexpiLdt", Kernel_Dimension::F);
    Matr       = S->reg<num_complex>("tmp.complex.Matr", Kernel_Dimension::FF);
}

int Kernel_Update_rho::exec_kernel_impl(int stat) {
    switch (Kernel_Representation::ele_repr_type) {
        case RepresentationPolicy::Diabatic: {
            for (int i = 0; i < Kernel_Dimension::F; ++i)
                invexpiEdt[i] = cos(E[i] * dt) - phys::math::im * sin(E[i] * dt);
            ARRAY_MATMUL3_TRANS2(U, T, invexpiEdt, T, Kernel_Dimension::F, Kernel_Dimension::F, 0, Kernel_Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, U, Kernel_Elec::rho_ele, U, Kernel_Dimension::F,
                                 Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, U, Kernel_Elec::rho_nuc, U, Kernel_Dimension::F,
                                 Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
            break;
        }
        case RepresentationPolicy::Adiabatic: {
            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
                ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
                ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
            }
            for (int i = 0; i < Kernel_Dimension::F; ++i)
                invexpiLdt[i] = cos(L[i] * dt) - phys::math::im * sin(L[i] * dt);
            ARRAY_MATMUL3_TRANS2(U, R, invexpiLdt, R, Kernel_Dimension::F, Kernel_Dimension::F, 0, Kernel_Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, U, Kernel_Elec::rho_ele, U, Kernel_Dimension::F,
                                 Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, U, Kernel_Elec::rho_nuc, U, Kernel_Dimension::F,
                                 Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F);
            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
                ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
                ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Kernel_Dimension::F, Kernel_Dimension::F, Kernel_Dimension::F,
                                     Kernel_Dimension::F);
            }
            break;
        }
        default:  // representation_policy::force, representation_policy::density
                  // LOG(FATAL);
            break;
    }
    return 0;
}

};  // namespace PROJECT_NS

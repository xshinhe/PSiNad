#include "Kernel_Update.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
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
    x              = DS->reg<num_real>("integrator.x", Dimension::N);
    p              = DS->reg<num_real>("integrator.p", Dimension::N);
    m              = DS->reg<num_real>("integrator.m", Dimension::N);
    minv           = DS->reg<num_real>("integrator.minv", Dimension::N);
    num_real* mass = DS->reg<num_real>("model.mass", Dimension::N);
    for (int i = 0; i < Dimension::N; ++i) { m[i] = mass[i], minv[i] = 1 / m[i]; }
}

int Kernel_Update_x::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::N; ++i) x[i] += p[i] * minv[i] * sdt;
    // ARRAY_SHOW(x, 1, Dimension::N);
    return 0;
}

void Kernel_Update_p::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);  //
    sdt = scale * dt;
}

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    f = DS->reg<num_real>("integrator.f", Dimension::N);
    p = DS->reg<num_real>("integrator.p", Dimension::N);
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::N; ++i) p[i] -= f[i] * sdt;
    // ARRAY_SHOW(p, 1, Dimension::N);
    return 0;
}



void Kernel_Update_T::read_param_impl(Param* PM) {
    dt     = PM->get<double>("dt", LOC(), phys::time_d);  //
    gammal = PM->get<double>("gammal", LOC(), 0.1);       //
    sdt    = scale * dt;
}

void Kernel_Update_T::init_data_impl(DataSet* S) {
    m = S->reg<num_real>("integrator.m", Dimension::N);
    p = S->reg<num_real>("integrator.p", Dimension::N);

    // if Langevin dynamics, set optimal c1 & c2p
    c1  = S->reg<num_real>("integrator.c1", Dimension::N);
    c2p = S->reg<num_real>("integrator.c2p", Dimension::N);
    for (int i = 0; i < Dimension::N; ++i) {
        c1[i]  = exp(-gammal * dt);
        c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
    }

    // if for NHC; registeration for auxiliary variables
    // ...
}

int Kernel_Update_T::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::N; ++i) {
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
    V          = S->reg<num_real>("model.V", Dimension::FF);
    E          = S->reg<num_real>("model.rep.E", Dimension::F);
    T          = S->reg<num_real>("model.rep.T", Dimension::FF);
    L          = S->reg<num_real>("model.rep.L", Dimension::F);
    R          = S->reg<num_complex>("model.rep.R", Dimension::FF);
    H          = S->reg<num_complex>("model.rep.H", Dimension::FF);
    U          = S->reg<num_complex>("integrator.U", Dimension::FF);
    invexpiEdt = S->reg<num_complex>("integrator.invexpiEdt", Dimension::F);
    invexpiLdt = S->reg<num_complex>("integrator.invexpiLdt", Dimension::F);
    Matr       = S->reg<num_complex>("tmp.complex.Matr", Dimension::FF);
}

int Kernel_Update_rho::exec_kernel_impl(int stat) {
    switch (Kernel_Representation::ele_repr_type) {
        case RepresentationPolicy::Diabatic: {
            for (int i = 0; i < Dimension::F; ++i) invexpiEdt[i] = cos(E[i] * dt) - phys::math::im * sin(E[i] * dt);
            ARRAY_MATMUL3_TRANS2(U, T, invexpiEdt, T, Dimension::F, Dimension::F, 0, Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, U, Kernel_Elec::rho_ele, U, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, U, Kernel_Elec::rho_nuc, U, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
            break;
        }
        case RepresentationPolicy::Adiabatic: {
            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
                ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T,  //
                                     Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                ARRAY_MATMUL3_TRANS1(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Dimension::F, Dimension::F, Dimension::F, Dimension::F);
            }
            for (int i = 0; i < Dimension::F; ++i) invexpiLdt[i] = cos(L[i] * dt) - phys::math::im * sin(L[i] * dt);
            ARRAY_MATMUL3_TRANS2(U, R, invexpiLdt, R, Dimension::F, Dimension::F, 0, Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, U, Kernel_Elec::rho_ele, U, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
            ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, U, Kernel_Elec::rho_nuc, U, Dimension::F, Dimension::F,
                                 Dimension::F, Dimension::F);
            if (Kernel_Representation::ini_repr_type == RepresentationPolicy::Diabatic) {
                ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_ele, T, Kernel_Elec::rho_ele, T,  //
                                     Dimension::F, Dimension::F, Dimension::F, Dimension::F);
                ARRAY_MATMUL3_TRANS2(Kernel_Elec::rho_nuc, T, Kernel_Elec::rho_nuc, T,  //
                                     Dimension::F, Dimension::F, Dimension::F, Dimension::F);
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

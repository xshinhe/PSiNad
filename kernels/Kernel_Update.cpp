#include "Kernel_Update.h"

#include "../core/linalg.h"
#include "Kernel_Declare.h"
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
    x              = DS->reg<num_real>("integrator.x", Dimension::PN);
    p              = DS->reg<num_real>("integrator.p", Dimension::PN);
    m              = DS->reg<num_real>("integrator.m", Dimension::PN);
    minv           = DS->reg<num_real>("integrator.minv", Dimension::PN);
    num_real* mass = DS->reg<num_real>("model.mass", Dimension::N);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_real* m    = this->m + iP * Dimension::N;
        num_real* minv = this->minv + iP * Dimension::N;
        for (int j = 0; j < Dimension::N; ++j) {
            m[j]    = mass[j];
            minv[j] = 1 / m[j];
        }
    }
}

int Kernel_Update_x::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::PN; ++i) x[i] += p[i] * minv[i] * sdt;
    return 0;
}

void Kernel_Update_p::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);
    sdt = scale * dt;
}

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    f = DS->reg<num_real>("integrator.f", Dimension::PN);
    p = DS->reg<num_real>("integrator.p", Dimension::PN);
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::PN; ++i) p[i] -= f[i] * sdt;
    return 0;
}

void Kernel_Update_T::read_param_impl(Param* PM) {
    dt     = PM->get<double>("dt", LOC(), phys::time_d);
    gammal = PM->get<double>("gammal", LOC(), 0.1);
    sdt    = scale * dt;
}

void Kernel_Update_T::init_data_impl(DataSet* S) {
    m = S->reg<num_real>("integrator.m", Dimension::PN);
    p = S->reg<num_real>("integrator.p", Dimension::PN);

    // if Langevin dynamics, set optimal c1 & c2p
    c1  = S->reg<num_real>("integrator.c1", Dimension::PN);
    c2p = S->reg<num_real>("integrator.c2p", Dimension::PN);
    for (int i = 0; i < Dimension::PN; ++i) {
        c1[i]  = exp(-gammal * dt);
        c2p[i] = sqrt(1.0 - c1[i] * c1[i]);
    }

    // if for NHC; registeration for auxiliary variables
    // ...
}

int Kernel_Update_T::exec_kernel_impl(int stat) {
    for (int i = 0; i < Dimension::PN; ++i) {
        Kernel_Random::rand_gaussian(&randu);
        p[i] = c1[i] * p[i] + c2p[i] * sqrt(m[i] / beta) * randu;
    }
    return 0;
}


void Kernel_Update_c::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);  //
    sdt = scale * dt;
}

void Kernel_Update_c::init_data_impl(DataSet* S) {
    E   = S->reg<num_real>("model.rep.E", Dimension::PF);
    T   = S->reg<num_real>("model.rep.T", Dimension::PFF);
    L   = S->reg<num_real>("model.rep.L", Dimension::PF);
    R   = S->reg<num_complex>("model.rep.R", Dimension::PFF);
    U   = S->reg<num_complex>("integrator.U", Dimension::PFF);
    Udt = S->reg<num_complex>("integrator.Udt", Dimension::PFF);

    invexpidiagdt = S->reg<num_complex>("integrator.tmp.invexpidiagdt", Dimension::F);
}

int Kernel_Update_c::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        // local variables for iP-th of swarm
        num_real* E      = this->E + iP * Dimension::F;
        num_real* T      = this->T + iP * Dimension::FF;
        num_real* L      = this->L + iP * Dimension::F;
        num_complex* R   = this->R + iP * Dimension::FF;
        num_complex* U   = this->U + iP * Dimension::FF;
        num_complex* Udt = this->Udt + iP * Dimension::FF;

        switch (Kernel_Representation::ele_repr_type) {
            case RepresentationPolicy::Diabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * E[i] * dt);
                ARRAY_MATMUL3_TRANS2(Udt, T, invexpidiagdt, T, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            case RepresentationPolicy::Adiabatic: {
                for (int i = 0; i < Dimension::F; ++i) invexpidiagdt[i] = exp(-phys::math::im * L[i] * dt);
                ARRAY_MATMUL3_TRANS2(Udt, R, invexpidiagdt, R, Dimension::F, Dimension::F, 0, Dimension::F);
                ARRAY_MATMUL(U, Udt, U, Dimension::F, Dimension::F, Dimension::F);
                break;
            }
            default:  // representation_policy::force, representation_policy::density
                      // LOG(FATAL);
                break;
        }
    }
    return 0;
}
};  // namespace PROJECT_NS

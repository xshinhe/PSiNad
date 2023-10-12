#include "Kernel_Update_p.h"

#include "Kernel_Declare.h"

namespace PROJECT_NS {

void Kernel_Update_p::read_param_impl(Param* PM) {
    dt  = PM->get<double>("dt", LOC(), phys::time_d);
    sdt = scale * dt;
}

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    f    = DS->reg<num_real>("integrator.f", Dimension::PN);
    p    = DS->reg<num_real>("integrator.p", Dimension::PN);
    minv = DS->reg<num_real>("integrator.minv", Dimension::PN);
    Ekin = DS->reg<num_real>("integrator.Ekin", Dimension::P);
}

void Kernel_Update_p::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_real* p    = this->p + iP * Dimension::N;
        num_real* minv = this->minv + iP * Dimension::N;
        num_real* Ekin = this->Ekin + iP;
        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
    }
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        num_real* f    = this->f + iP * Dimension::N;
        num_real* p    = this->p + iP * Dimension::N;
        num_real* minv = this->minv + iP * Dimension::N;
        num_real* Ekin = this->Ekin + iP;

        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) {
            p[i] -= f[i] * sdt;
            Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
        }
    }
    return 0;
}

};  // namespace PROJECT_NS
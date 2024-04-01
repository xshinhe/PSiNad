#include "Kernel_Update_p.h"

#include "../core/vars_list.h"
#include "Kernel_Declare.h"

namespace PROJECT_NS {

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    dt_ptr   = DS->def(DATA::iter::dt);
    f        = DS->def(DATA::integrator::f);
    p        = DS->def(DATA::integrator::p);
    minv     = DS->def(DATA::integrator::minv);
    Ekin     = DS->def(DATA::integrator::Ekin);
    frez_ptr = DS->def(DATA::iter::frez);
}

void Kernel_Update_p::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;
        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
    }
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    if (frez_ptr[0]) return 0;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* f    = this->f + iP * Dimension::N;
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;

        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) {
            p[i] -= f[i] * scale * dt_ptr[0];
            Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
        }
    }
    return 0;
}

};  // namespace PROJECT_NS
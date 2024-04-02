#include "Kernel_Update_x.h"

#include "../core/vars_list.h"

namespace PROJECT_NS {

void Kernel_Update_x::init_data_impl(DataSet* DS) {
    dt_ptr          = DS->def(DATA::iter::dt);
    x               = DS->def(DATA::integrator::x);
    p               = DS->def(DATA::integrator::p);
    m               = DS->def(DATA::integrator::m);
    minv            = DS->def(DATA::integrator::minv);
    kids_real* mass = DS->def(DATA::model::mass);
    frez_ptr        = DS->def(DATA::iter::frez);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* m    = this->m + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        for (int j = 0; j < Dimension::N; ++j) {
            m[j]    = mass[j];
            minv[j] = 1 / m[j];
        }
    }
}

int Kernel_Update_x::exec_kernel_impl(int stat) {
    if (frez_ptr[0]) return 0;
    for (int i = 0; i < Dimension::PN; ++i) x[i] += p[i] * minv[i] * scale * dt_ptr[0];
    return 0;
}

};  // namespace PROJECT_NS

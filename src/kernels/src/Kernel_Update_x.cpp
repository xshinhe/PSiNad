#include "kids/Kernel_Update_x.h"

#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_x::getName() { return "Kernel_Update_x"; }

int Kernel_Update_x::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_x::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    x               = DS->def(DATA::integrator::x);
    p               = DS->def(DATA::integrator::p);
    m               = DS->def(DATA::integrator::m);
    minv            = DS->def(DATA::integrator::minv);
    kids_real* mass = DS->def(DATA::model::mass);
    dt              = DS->def(DATA::flowcontrol::dt);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* m    = this->m + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        for (int j = 0; j < Dimension::N; ++j) {
            m[j]    = mass[j];
            minv[j] = 1 / m[j];
        }
    }
}

Status& Kernel_Update_x::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;
    for (int i = 0; i < Dimension::PN; ++i) x[i] += p[i] * minv[i] * scale * dt[0];
    return stat;
}

};  // namespace PROJECT_NS

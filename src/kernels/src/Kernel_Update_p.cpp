#include "kids/Kernel_Update_p.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_p::getName() { return "Kernel_Update_p"; }

int Kernel_Update_p::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_p::setInputDataSet_impl(std::shared_ptr<DataSet>& DS) {
    dt_ptr   = DS->def(DATA::iter::dt);
    f        = DS->def(DATA::integrator::f);
    p        = DS->def(DATA::integrator::p);
    minv     = DS->def(DATA::integrator::minv);
    Ekin     = DS->def(DATA::integrator::Ekin);
    frez_ptr = DS->def(DATA::iter::frez);
}

Status& Kernel_Update_p::initializeKernel_impl(Status& stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;
        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
    }
    return stat;
}

Status& Kernel_Update_p::executeKernel_impl(Status& stat) {
    if (frez_ptr[0]) return stat;

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
    return stat;
}

};  // namespace PROJECT_NS
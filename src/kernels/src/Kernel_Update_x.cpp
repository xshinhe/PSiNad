#include "kids/Kernel_Update_x.h"

#include "kids/Kernel_Monodromy.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Update_x::getName() { return "Kernel_Update_x"; }

int Kernel_Update_x::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Update_x::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    x         = DS->def(DATA::integrator::x);
    p         = DS->def(DATA::integrator::p);
    m         = DS->def(DATA::integrator::m);
    minv      = DS->def(DATA::integrator::minv);
    mono      = DS->def(DATA::integrator::monodromy::mono);
    monodt    = DS->def(DATA::integrator::monodromy::monodt);
    auto mass = DS->def(DATA::model::mass);
    dt        = DS->def(DATA::flowcontrol::dt);
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto m    = this->m.subspan(iP * Dimension::N, Dimension::N);
        auto minv = this->minv.subspan(iP * Dimension::N, Dimension::N);
        for (int j = 0; j < Dimension::N; ++j) {
            m[j]    = mass[j];
            minv[j] = 1 / m[j];
        }
    }
}

Status& Kernel_Update_x::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;
    for (int i = 0; i < Dimension::PN; ++i) x[i] += p[i] * minv[i] * scale * dt[0];
    // trace on monodromy
    if (Kernel_Monodromy::enable) update_monodromy();
    return stat;
}

void Kernel_Update_x::update_monodromy() {
    int N0   = 0;
    int N1   = Dimension::N;
    int N2   = Dimension::N + Dimension::F;
    int N3   = 2 * Dimension::N + Dimension::F;
    int N4   = 2 * Dimension::N + 2 * Dimension::F;
    int N4N4 = N4 * N4;
    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto mono   = this->mono.subspan(iP * N4N4, N4N4);    //
        auto monodt = this->monodt.subspan(iP * N4N4, N4N4);  //
        ARRAY_EYE(monodt.data(), N4);
        for (int I = 0; I < Dimension::N; ++I) {
            for (int J = 0; J < Dimension::N; ++J) { monodt[(N0 + I) * N4 + (N2 + J)] = minv[J] * scale * dt[0]; }
        }
        ARRAY_MATMUL(mono.data(), monodt.data(), mono.data(), N4, N4, N4);
    }
}

};  // namespace PROJECT_NS

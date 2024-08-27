#include "kids/Kernel_Conserve.h"

#include <algorithm>

#include "kids/Kernel_Representation.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Conserve::getName() { return "Kernel_Conserve"; }

int Kernel_Conserve::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Conserve::setInputParam_impl(std::shared_ptr<Param> PM) {
    conserve_scale = _param->get_bool({"solver.conserve_scale"}, LOC(), false);  // default as false
}

void Kernel_Conserve::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    vpes             = DS->def(DATA::model::vpes);
    Ekin             = DS->def(DATA::integrator::Ekin);
    Epot             = DS->def(DATA::integrator::Epot);
    Etot             = DS->def(DATA::integrator::Etot);
    Etot_init        = DS->def(DATA::init::Etot);
    Etot_prev        = DS->def(DATA::last::Etot);
    E                = DS->def(DATA::model::rep::E);
    m                = DS->def(DATA::integrator::m);
    p                = DS->def(DATA::integrator::p);
    succ_ptr         = DS->def(DATA::flowcontrol::succ);
    last_attempt_ptr = DS->def(DATA::flowcontrol::last_attempt);
    frez_ptr         = DS->def(DATA::flowcontrol::frez);
    fail_type_ptr    = DS->def(DATA::flowcontrol::fail_type);
}

Status& Kernel_Conserve::initializeKernel_impl(Status& stat) {
    bool conserve_scale_bak = conserve_scale;
    conserve_scale          = false;
    // executeKernel(stat);  // calc Epot
    conserve_scale = conserve_scale_bak;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        Etot[iP] = Ekin[iP] + Epot[iP], Etot_init[iP] = Etot[iP];
        Etot_prev[iP] = Etot_init[iP];
    }
    return stat;
};

#define DECLARE_LOCAL_SPAN(varname, iparrallel, size) auto varname = this->varname.subspan(iparrallel * size, size)

Status& Kernel_Conserve::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        // DECLARE_LOCAL_SPAN(E, iP, Dimension::F);

        auto E         = this->E.subspan(iP * Dimension::F, Dimension::F);
        auto p         = this->p.subspan(iP * Dimension::N, Dimension::N);
        auto m         = this->m.subspan(iP * Dimension::N, Dimension::N);
        auto Etot      = this->Etot.subspan(iP, 1);
        auto Etot_init = this->Etot_init.subspan(iP, 1);
        auto Etot_prev = this->Etot_prev.subspan(iP, 1);
        auto Ekin      = this->Ekin.subspan(iP, 1);
        auto Epot      = this->Epot.subspan(iP, 1);
        auto vpes      = this->vpes.subspan(iP, 1);
        // kids_real* E         = this->E + iP * Dimension::F;
        // kids_real* p         = this->p + iP * Dimension::N;
        // kids_real* m         = this->m + iP * Dimension::N;
        // kids_real* Etot      = this->Etot + iP;
        // kids_real* Etot_init = this->Etot_init + iP;
        // kids_real* Etot_prev = this->Etot_prev + iP;
        // kids_real* Ekin      = this->Ekin + iP;
        // kids_real* Epot      = this->Epot + iP;
        // kids_real* vpes      = this->vpes + iP;

        Ekin[0] = 0.0e0;
        for (int j = 0; j < Dimension::N; ++j) Ekin[0] += 0.5e0 * p[j] * p[j] / m[j];

        if (Kernel_Representation::onthefly) {
            double thres = 0.02;
            if (stat.last_attempt && stat.fail_type != 2) { cnt_loose = 3; }
            if (cnt_loose > 0) {
                thres = cnt_loose * cnt_loose;  // help for last_attempt of mndo failure
                std::cout << "try loose thres = " << thres << " because failure of QM\n";
                cnt_loose--;
            }

            if (stat.last_attempt && stat.fail_type == 2) {
                thres = 0.5;  // loose threshold
                std::cout << "try loose thres = " << thres << " because of failure of CONSERVE\n";
            }

            if (std::abs(Ekin[0] + Epot[0] - Etot_prev[0]) * phys::au_2_kcal_1mea > thres) {
                std::cout << "fail in conserve ERROR: "  //
                          << fabs(Ekin[0] + Epot[0] - Etot_prev[0]) * phys::au_2_kcal_1mea << " > " << thres << "\n";
                stat.succ      = false;
                stat.fail_type = 2;
            } else {
                Etot_prev[0]   = Ekin[0] + Epot[0];
                stat.succ      = true;
                stat.fail_type = 0;
            }
        }

        if (conserve_scale) {
            double scale = std::sqrt(std::max({Etot_init[0] - Epot[0], 0.0e0}) / Ekin[0]);
            for (int j = 0; j < Dimension::N; ++j) p[j] *= scale;
            // Ekin[0] = Etot_init[0] - Epot[0];
        }
        Etot[0] = Epot[0] + Ekin[0];
    }
    return stat;
}


};  // namespace PROJECT_NS

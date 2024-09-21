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
    conserve_scale   = _param->get_bool({"solver.conserve_scale"}, LOC(), false);  // default as false
    thres_kcalpermol = _param->get_real({"solver.thres_kcalpermol"}, LOC(), 0.02);
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

Status& Kernel_Conserve::executeKernel_impl(Status& stat) {
    if (stat.frozen) return stat;
    if (!stat.succ) return stat;  // reject conservation for un-succesful calculation

    for (int iP = 0; iP < Dimension::P; ++iP) {
        auto E         = this->E.subspan(iP * Dimension::F, Dimension::F);
        auto p         = this->p.subspan(iP * Dimension::N, Dimension::N);
        auto m         = this->m.subspan(iP * Dimension::N, Dimension::N);
        auto Etot      = this->Etot.subspan(iP, 1);
        auto Etot_init = this->Etot_init.subspan(iP, 1);
        auto Etot_prev = this->Etot_prev.subspan(iP, 1);
        auto Ekin      = this->Ekin.subspan(iP, 1);
        auto Epot      = this->Epot.subspan(iP, 1);
        auto vpes      = this->vpes.subspan(iP, 1);

        Ekin[0] = 0.0e0;
        for (int j = 0; j < Dimension::N; ++j) Ekin[0] += 0.5e0 * p[j] * p[j] / m[j];

        if (Kernel_Representation::onthefly) {
            double thres = thres_kcalpermol;

            // if try last QM calculation. we totally loose energy conservation in it's first several steps
            const int nstep_loose = 5;
            if (stat.last_attempt && stat.fail_type == 1) { cnt_loose = nstep_loose; }
            if (stat.fail_type != 2 && cnt_loose > 0) {
                thres = (cnt_loose * 5.0e0 + (nstep_loose - cnt_loose) * thres_kcalpermol) / nstep_loose;
                std::cout << "try loose thres = " << thres << " because failure of QM\n";
                cnt_loose--;
            }

            bool loose_10 = false;
            if (stat.last_attempt && stat.fail_type == 2) {
                loose_10 = true;
                thres = 10.0 * thres_kcalpermol;  // loose threshold
                std::cout << "try loose thres = " << thres << " because of failure of CONSERVATION\n";
            }

            double deltaE = fabs(Ekin[0] + Epot[0] - Etot_prev[0]) * phys::au_2_kcal_1mea;
            if (deltaE > thres) {
                std::cout << "fail in conserve ERROR: "                                     //
                          << deltaE  //
                          << " > " << thres << "\n";
                stat.succ      = false;
                stat.fail_type = 2;
                if(stat.first_step){
                    // the first step fail in energy conversation, so kinematic energy to too large for dt suggested
                    std::cout << "warning: the conservation for the first step fails! but we continue to run with loose of threshold\n";
                    std::cout << "but you'd better kill job and check the initial condition\n";
                    thres_kcalpermol = deltaE * 1.01;
                }
                if(loose_10 && Etot_prev[0] > Epot[0]) {
                    std::cout << "force scale the energy and recover the trajectory! it should be carefull!!!" << std::endl;
                    double scale = std::sqrt(std::max({Etot_prev[0] - Epot[0], 0.0e0}) / Ekin[0]);
                    for (int j = 0; j < Dimension::N; ++j) p[j] *= scale;
                    stat.succ      = true;
                    stat.fail_type = 0;
                }
            } else {
                std::cout << "now deltaE: "                                     //
                          << deltaE  //
                          << " <= " << thres << "\n";
                Etot_prev[0]   = Ekin[0] + Epot[0];
                stat.fail_type = 0;  // as long as both succeed, we reset the fail_type as 0
                //if(!stat.last_attempt && thres / deltaE < 1.414e0) { // check if it is reset to -1
                //    stat.fail_type = -1; // succ but dangerous, so don't increase the step size, just keep it
                //}
            }
        }

        if (conserve_scale) {
            double scale = std::sqrt(std::max({Etot_init[0] - Epot[0], 0.0e0}) / Ekin[0]);
            for (int j = 0; j < Dimension::N; ++j) p[j] *= scale;
        }
        Etot[0] = Epot[0] + Ekin[0];
    }
    return stat;
}


};  // namespace PROJECT_NS

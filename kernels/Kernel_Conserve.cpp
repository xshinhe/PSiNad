#include "Kernel_Conserve.h"

#include "Kernel_Declare.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

void Kernel_Conserve::read_param_impl(Param* PM) {
    conserve_scale = PM->get<bool>("conserve_scale", LOC(), false);  // default as false
}

void Kernel_Conserve::init_data_impl(DataSet* DS) {
    vpes      = DS->def<kids_real>("model.vpes", Dimension::P);
    Ekin      = DS->def<kids_real>("integrator.Ekin", Dimension::P);
    Epot      = DS->def<kids_real>("integrator.Epot", Dimension::P);
    Etot      = DS->def<kids_real>("integrator.Etot", Dimension::P);
    Etot_prev = DS->def<kids_real>("integrator.Etot_prev", Dimension::P);
    E         = DS->def<kids_real>("model.rep.E", Dimension::PF);
    m         = DS->def<kids_real>("integrator.m", Dimension::PN);
    p         = DS->def<kids_real>("integrator.p", Dimension::PN);
    succ_ptr  = DS->def<bool>("iter.succ");
}

void Kernel_Conserve::init_calc_impl(int stat) {
    Etot_init = _DataSet->def("init.Etot", Etot, Dimension::P);

    bool conserve_scale_bak = conserve_scale;
    conserve_scale          = false;
    exec_kernel(stat);  // calc Epot
    conserve_scale = conserve_scale_bak;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        Etot[iP] = Ekin[iP] + Epot[iP], Etot_init[iP] = Etot[iP];
        Etot_prev[iP] = Etot_init[iP];
    }
};

int Kernel_Conserve::exec_kernel_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* E         = this->E + iP * Dimension::F;
        kids_real* p         = this->p + iP * Dimension::N;
        kids_real* m         = this->m + iP * Dimension::N;
        kids_real* Etot      = this->Etot + iP;
        kids_real* Etot_init = this->Etot_init + iP;
        kids_real* Ekin      = this->Ekin + iP;
        kids_real* Epot      = this->Epot + iP;
        kids_real* vpes      = this->vpes + iP;

        // Epot[0] = vpes[0] + E[(*Kernel_Elec::occ_nuc)];
        Ekin[0] = 0.0e0;
        for (int j = 0; j < Dimension::N; ++j) Ekin[0] += 0.5e0 * p[j] * p[j] / m[j];

        // double normxx = 0.0;
        // for (int i = 0; i < Dimension::N; ++i) normxx += p[i] * p[i];
        // normxx = sqrt(normxx);
        // std::cout << "norm = " << normxx << "\n";
        // std::cout << "Ekin = " << Ekin[0] << "\n";

        bool prev_succ = succ_ptr[0];

        if (fabs(Ekin[0] + Epot[0] - Etot_prev[0]) * phys::au_2_kcal_1mea > 0.05) {
            std::cout << "ABS ERROR: " << fabs(Ekin[0] + Epot[0] - Etot_prev[0]) * phys::au_2_kcal_1mea << "\n";
            std::cout << "REL ERROR: " << (Ekin[0] + Epot[0]) / (Etot_prev[0]) << "\n";
            succ_ptr[0] = false;
        } else {
            Etot_prev[0] = Ekin[0] + Epot[0];
        }

        // std::cout << "Etot = " << Ekin[0] + Epot[0] << "\n";
        // std::cout << "2 % = " << (Ekin[0] + Epot[0]) / (Etot_prev[0]) << "\n";
        // std::cout << phys::au_2_kcal_1mea << "\n";

        if (prev_succ && conserve_scale) {
            double scale = std::sqrt(std::max({Etot_init[0] - Epot[0], 0.0e0}) / Ekin[0]);
            for (int j = 0; j < Dimension::N; ++j) p[j] *= scale;
            Ekin[0] = Etot_init[0] - Epot[0];
        }
    }
    return 0;
}


};  // namespace PROJECT_NS

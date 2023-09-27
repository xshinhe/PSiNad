#include "./interf_sander.h"

// used for qspcfw model temporarily

sander_ForceField::sander_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    CHECK_EQ(N % 9, 0);
    natom = N / 3;
    nmol  = N / 9;
    inpcrd = Param_GetT(std::string, parm, "inpcrd", "inpcrd");
    prmtop = Param_GetT(std::string, parm, "prmtop", "prmtop");

    const double m_H = ELEMENTS_MASS_NOAVG[1] / phys::au_2_amu, m_C = ELEMENTS_MASS_NOAVG[6] / phys::au_2_amu,
          m_N = ELEMENTS_MASS_NOAVG[7] / phys::au_2_amu, m_O = ELEMENTS_MASS_NOAVG[8] / phys::au_2_amu;

    const double mass_h2o[9] = {m_O, m_O, m_O, m_H, m_H, m_H, m_H, m_H, m_H};

    int ierr;
    double box[6];
    ierr = read_inpcrd_file(inpcrd.c_str(), mod_R0, box);
    pme_sander_input(&inp);
    inp.saltcon = 0.2;
    inp.cut = Param_GetT(double, parm, "cut", 7.0);
    ierr = sander_setup_mm(prmtop.c_str(), mod_R0, box, &inp);


    crd = new double[N];
    frc = new double[N];

    for (int i = 0; i < nmol; ++i)
        for (int j = 0; j < 9; ++j) {
            mod_M[i*9 + j] = mass_h2o[j];
            mod_R0[i*9 + j] /= phys::au_2_ang;
        }

    tag = name() + tag;
}

sander_ForceField::~sander_ForceField() {
    sander_cleanup();
    delete [] crd;
    delete [] frc;
};

int sander_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    return 0;
}

int sander_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
        const int& rdim) {
    for (int i = 0; i < N; ++i) crd[i] = R[i] * phys::au_2_ang;
    set_positions(crd);
    energy_forces(&ene, frc);
    *V = ene.tot / phys::au_2_kcal_1mea;
    //*V = ene.tot;

    for (int i = 0; i < N; ++i) {
        dV[i] = -frc[i] * phys::au_2_ang / phys::au_2_kcal_1mea;
    }

    return 0;
}



#include "smallmol_models.h"

#include "../../utils/elements.h"
/**
 * @brief interface to fortran H2O forcefield
 *
 */
extern "C" {
void h2opot_ccall(double*, double*, double*, double*, int*);
};
int ForceField_npes_H2O(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                        const int& rdim) {  // @todo add box?
    int flagcopy = flag;
    h2opot_ccall(V, dV, ddV, R, &flagcopy);
    return 0;
}

/**
 * @brief interface to fortran H2O2 forcefield
 *
 * @param x position array (size=12)
 * @param dv force array (size=12)
 * @param v energy
 * @param vb magnitude of dv, should always be 1.0
 * @todo clearify above parameters
 */
extern "C" {  // from fortran function pes_h2o2_mp_surf_b_
void pes_h2o2_surf_b(double* x, double* dv, double* v, double* vb);
void pes_h2o2_surf(double* x, double* v);
};
int ForceField_npes_H2O2(double* V, double* dV, double* ddV, double* R, double* P, const int& flag, const int& rdim) {
    double vb = 1.0;
    for (int i = 0; i < 12; ++i) R[i] *= phys::au_2_ang;
    pes_h2o2_surf_b(R, dV, V, &vb);
    for (int i = 0; i < 12; ++i) R[i] /= phys::au_2_ang;
    for (int i = 0; i < 12; ++i) dV[i] *= phys::au_2_ang;
    return 0;
}

/**
 * @brief interface to fortran NH3 forcefield
 *
 */
extern "C" {
void nh3_pot_cart_ccall(double*, double*, double*, double*, int*);
};
int ForceField_npes_NH3(double* V, double* dV, double* ddV, double* R, double* P, const int& flag, const int& rdim) {
    int flagcopy = flag;
    nh3_pot_cart_ccall(V, dV, ddV, R, &flagcopy);
    return 0;
}


/**
 * @brief interface to fortran CH2O forcefield
 *
 */
extern "C" {};
int ForceField_npes_CH2O(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                         const int& rdim) {  // @todo add box?
    // needed codes here
    return 0;
}

/**
 * @brief interface to para-H2 forcefield
 *
 */
int paraH2_pot(double* V, double* dV, double* nr, const int& N, const double& BoxL);
int ForceField_npes_paraH2(double* V, double* dV, double* ddV, double* R, double* P, const int& flag, const int& rdim,
                           const double& L) {
    paraH2_pot(V, dV, R, rdim, L);
    return 0;
}


SmallMol_ForceField::SmallMol_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    CHECK_EQ(N % 3, 0);
    Natom = N / 3;

    std::string ffflag = Param_GetT(std::string, parm, "type");
    try {
        fftype = SmallMol_Policy::_dict.at(ffflag);
    } catch (std::out_of_range& e) { LOG(FATAL) << "Unknown /type/ " << ffflag; }

    const double m_H = ELEMENTS_MASS_NOAVG[1] / phys::au_2_amu, m_C = ELEMENTS_MASS_NOAVG[6] / phys::au_2_amu,
                 m_N = ELEMENTS_MASS_NOAVG[7] / phys::au_2_amu, m_O = ELEMENTS_MASS_NOAVG[8] / phys::au_2_amu;
    switch (fftype) {
        case SmallMol_Policy::SM_H2O: {
            CHECK_EQ(N, 9);
            Nmole                    = 1;
            const double H2O_mass[9] = {m_H, m_H, m_H, m_O, m_O, m_O, m_H, m_H, m_H};
            const double H2O_eq[9]   = {1.77, 0.39, -0.525, 0.00646, 0.067, -0.016, -0.606, 1.68, 0.44};
            for (int i = 0; i < N; ++i) mod_R0[i] = H2O_eq[i], mod_M[i] = H2O_mass[i];
            break;
        }
        case SmallMol_Policy::SM_H2O2: {
            CHECK_EQ(N, 12);
            Nmole                      = 1;
            const double H2O2_mass[12] = {m_H, m_H, m_H, m_O, m_O, m_O, m_O, m_O, m_O, m_H, m_H, m_H};
            const double H2O2_eq[12]   = {0.94675,   -0.0589902, 0.946786,  0.0210381,  0.0912028, 0.724262,
                                        0.0912028, 0.0210381,  -0.724262, -0.0589902, 0.94675,   -0.946786};
            for (int i = 0; i < N; ++i) mod_R0[i] = H2O2_eq[i] / phys::au_2_ang, mod_M[i] = H2O2_mass[i];
            break;
        }
        case SmallMol_Policy::SM_NH3: {
            CHECK_EQ(N, 12);
            Nmole                     = 1;
            const double NH3_mass[12] = {// to be check later
                                         m_N, m_N, m_N, m_H, m_H, m_H, m_H, m_H, m_H, m_H, m_H, m_H};
            const double NH3_eq[12]   = {};
            for (int i = 0; i < N; ++i) mod_R0[i] = NH3_eq[i], mod_M[i] = NH3_mass[i];
            break;
        }
        case SmallMol_Policy::SM_CH2O: {
            CHECK_EQ(N, 12);
            Nmole = 1;
            // mass need to be
            const double CH2O_eq[12] = {};
            for (int i = 0; i < N; ++i) mod_R0[i] = CH2O_eq[i];
            break;
        }
        case SmallMol_Policy::SM_paraH2: {
            // CHECK_EQ(Natom % 2, 0);
            Nmole       = Natom;
            boxL        = Param_GetT(double, parm, "boxL", -1.0f);  // get box size (required!!!)
            double dens = Param_GetQ(phys::mass_density_d, parm, "dens", -1.0f);
            if (boxL < 0.0f) { Param_Reset(boxL, pow(Natom * m_H * 2 / dens, 1.0f / 3.0f)); }
            // std::cout << boxL <<std::endl;
            CHECK_GT(boxL, 0);
            for (int i = 0; i < N; ++i) mod_M[i] = m_H * 2;

            int N_box  = int(pow(Natom * 1.0f, 1.0f / 3.0f)) + 1;
            double L_r = boxL / N_box;
            int a1, a2, a3;
            int idxR = 0;
            for (int iatom = 1; iatom < Natom + 1; ++iatom) {
                a1             = floor(iatom / (N_box * N_box));
                a2             = floor((iatom - a1 * (N_box * N_box)) / N_box);
                a3             = iatom % N_box;
                mod_R0[idxR++] = a1 * L_r;
                mod_R0[idxR++] = a2 * L_r;
                mod_R0[idxR++] = a3 * L_r;
                // std::cout << a1 << a2 << a3 <<std::endl;
            }
            break;
        }
        default:
            LOG(FATAL);
    }
    for (int i = 0; i < N; ++i) { mod_P0[i] = 0.0f, mod_sigmaR[i] = 0.0f, mod_sigmaP[i] = 0.0f; }

    tag = name() + ffflag + "_" + tag;
}

SmallMol_ForceField::~SmallMol_ForceField(){};

// use general init function
int SmallMol_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    return 0;
}

// use default specificatin function
int SmallMol_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return BO_ForceField ::ForceField_spec(nr, np, nm, rdim);
}

// call pes
int SmallMol_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                         const int& rdim) {
    switch (fftype) {
        case SmallMol_Policy::SM_H2O:
            ForceField_npes_H2O(V, dV, ddV, R, P, flag, rdim);
            break;
        case SmallMol_Policy::SM_H2O2:
            ForceField_npes_H2O2(V, dV, ddV, R, P, flag, rdim);
            break;
        case SmallMol_Policy::SM_NH3:
            ForceField_npes_NH3(V, dV, ddV, R, P, flag, rdim);
            break;
        case SmallMol_Policy::SM_CH2O:
            ForceField_npes_CH2O(V, dV, ddV, R, P, flag, rdim);
            break;
        case SmallMol_Policy::SM_paraH2:
            ForceField_npes_paraH2(V, dV, ddV, R, P, flag, rdim, boxL);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

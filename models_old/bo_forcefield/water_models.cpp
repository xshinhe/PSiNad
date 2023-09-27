#include "water_models.h"

#include "../../utils/elements.h"

/**
 * @description:fortran interface of qspcfw water model
 * @param Natom number of atoms
 * @param pbox box length array(3)
 * @param coord coordinate of {O, H, H}
 * @param charge charge of each atoms
 * @param force output force
 * @param energy output potential energy
 * @param cutoff cutoff length, should be less than half of pbox
 * @param ewald_parm ewald summation, between 0 to 1
 * @return {*}
 */
extern "C" {
void qspcfw_force(int* Natom, double* pbox, double* coord, double* charge, double* force, double* energy,
                  double* cutoff, double* ewald_parm);
}

// const double au_2_kcal = 627.5094;
int Water_qspcfw_npes(double* V, double* dV, double* R, int& Natom, double* pbox, double* charge, double& cutoff,
                      double& ewald_parm) {
    for (int i = 0; i < Natom * 3; ++i) R[i] *= phys::au_2_ang;
    LOG(INFO) << dV[0];
    qspcfw_force(&Natom, pbox, R, charge, dV, V, &cutoff, &ewald_parm);
    // LOG(INFO) <<"ang per au"<< FMT(8) << phys::au_2_kcal_1mea;
    for (int i = 0; i < Natom * 3; ++i) {
        R[i] /= phys::au_2_ang;
        // std::cout << FMT(8) << dV[i];
        dV[i] *= -phys::au_2_ang / phys::au_2_kcal_1mea;
    }
    *V /= phys::au_2_kcal_1mea;
    return 0;
}

Water_ForceField::Water_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    CHECK_EQ(N % 3, 0);
    Natom = N / 3;
    CHECK_EQ(Natom % 3, 0);

    std::string ffflag = Param_GetT(std::string, parm, "type");
    try {
        fftype = WaterModelPolicy::_dict.at(ffflag);
    } catch (std::out_of_range& e) { LOG(FATAL) << "Unknown /type/ " << ffflag; }
    const double m_H = ELEMENTS_MASS_NOAVG[1] / phys::au_2_amu, m_O = ELEMENTS_MASS_NOAVG[8] / phys::au_2_amu;
    switch (fftype) {
        case WaterModelPolicy::qspcfw: {
            tag = "qspcfw";
            // get box size and cutoff from param.json
            boxL        = Param_GetQ(phys::length_d, parm, "boxL",
                              -1.0f);  // get box size (required!!!)
            cutoff      = Param_GetQ(phys::length_d, parm, "cutoff",
                                -1.0f);  // get cutoff length (required!!!)
            double dens = Param_GetQ(phys::mass_density_d, parm, "dens", -1.0f);
            if (boxL < 0.0f) Param_Reset(boxL, pow(Natom * m_H * 2 / dens, 1.0f / 3.0f));
            if (cutoff < 0.0f) {
                Param_Reset(cutoff, boxL / 2);
                LOG(INFO) << "cutoff distance reset to half of box length.";
            }

            CHECK_LE(cutoff, boxL / 2);

            ALLOCATE_PTR_TO_VECTOR(pbox, 3);
            pbox[0] = boxL * phys::au_2_ang;
            pbox[1] = boxL * phys::au_2_ang;
            pbox[2] = boxL * phys::au_2_ang;
            cutoff *= phys::au_2_ang;
            // get charge array
            double chargeO  = Param_GetQ(phys::electric_charge_d, parm, "chargeO", -1.0f);
            double chargeH1 = Param_GetQ(phys::electric_charge_d, parm, "chargeH1", 0.5f);
            double chargeH2 = Param_GetQ(phys::electric_charge_d, parm, "chargeH2", 0.5f);

            ALLOCATE_PTR_TO_VECTOR(charge_arr, Natom);
            for (int i = 0; i < Natom; i += 3) {
                charge_arr[i]     = chargeO;
                charge_arr[i + 1] = chargeH1;
                charge_arr[i + 2] = chargeH2;
            }
            // get ewald paramter
            ewald_parm = Param_GetT(double, parm, "ewald parameter", 0.0f);
            CHECK_GE(ewald_parm, 0.0f);
            CHECK_LE(ewald_parm, 1.0f);
            break;
        }
        default:
            LOG(FATAL);
            break;
    }
    for (int i = 0; i < N; ++i) { mod_P0[i] = 0.0f, mod_sigmaR[i] = 0.0f, mod_sigmaP[i] = 0.0f; }

    tag = name() + ffflag + "_" + tag;
}

int Water_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    return 0;
}

int Water_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return BO_ForceField ::ForceField_spec(nr, np, nm, rdim);
}

int Water_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                      const int& rdim) {
    return Water_qspcfw_npes(V, dV, R, Natom, pbox, charge_arr, cutoff, ewald_parm);
}

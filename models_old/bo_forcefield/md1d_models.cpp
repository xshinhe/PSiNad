#include "md1d_models.h"

#include "../../utils/definitions.h"


MD1D_ForceField::MD1D_ForceField(const Param& iparm, const int& child) : BO_ForceField(iparm){};  // for child to call
MD1D_ForceField::MD1D_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    LOG(WARNING) << iparm;

    CHECK_EQ(N, 1);
    double r0     = Param_GetV(r0, parm, 1.0f);
    double p0     = Param_GetV(p0, parm, 0.0f);
    mod_R0[0]     = r0;
    mod_P0[0]     = p0;
    mod_sigmaR[0] = 0.0f;
    mod_sigmaP[0] = 0.0f;

    ffflag = Param_GetT(std::string, parm, "ffflag", "#ho");
    fftype = MD1DPolicy::_dict.at(ffflag);

    switch (fftype) {
        case MD1DPolicy::MD1D_HO:
            mod_M[0] = 1.0f;
            mod_W[0] = 1.0f;
            break;
        default:
            LOG(FATAL) << "unsupported name";
    }

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

MD1D_ForceField::~MD1D_ForceField(){};

int MD1D_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    return 0;
}

int MD1D_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return BO_ForceField ::ForceField_spec(nr, np, nm, rdim);
}

int MD1D_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                     const int& rdim) {
    V[0] = 0.5 * mod_M[0] * mod_W[0] * mod_W[0] * R[0] * R[0];
    if (flag < 1) return 0;
    dV[0] = mod_M[0] * mod_W[0] * mod_W[0] * R[0];
    if (flag < 2) return 0;
    ddV[0] = mod_M[0] * mod_W[0] * mod_W[0];
    return 0;
}

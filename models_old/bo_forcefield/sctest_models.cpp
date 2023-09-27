#include "sctest_models.h"

#include "../../utils/definitions.h"

SCTEST_ForceField::SCTEST_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    LOG(WARNING) << iparm;

    double r0     = Param_GetV(r0, parm, 1.0f);
    double p0     = Param_GetV(p0, parm, 0.0f);
    mod_R0[0]     = r0;
    mod_P0[0]     = p0;
    mod_sigmaR[0] = 0.0f;
    mod_sigmaP[0] = 0.0f;

    ffflag = Param_GetT(std::string, parm, "ffflag", "#sc1d");
    fftype = SCTESTPolicy::_dict.at(ffflag);

    switch (fftype) {
        case SCTESTPolicy::SC1D:
            CHECK_EQ(N, 1);
            mod_M[0] = 1.0f;
            mod_W[0] = std::sqrt(2.0f);
            break;
        case SCTESTPolicy::SC2D:
            CHECK_EQ(N, 2);
            mod_M[0] = 1.0f;
            mod_W[0] = std::sqrt(2.0f);
            mod_M[1] = 1.0f;
            mod_W[1] = 1.0f / 3.0f;
            break;
        default:
            LOG(FATAL) << "unsupported name";
    }

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

SCTEST_ForceField::~SCTEST_ForceField(){};

int SCTEST_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& icycle) {
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, icycle);
    return 0;
}

int SCTEST_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return BO_ForceField ::ForceField_spec(nr, np, nm, rdim);
}


int SCTEST_ForceField::ForceField_npes_SC1D(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                            const int& rdim) {
    const double k      = mod_M[0] * mod_W[0] * mod_M[0];
    const double parm_a = -0.1f, parm_b = 0.1f;
    V[0]   = R[0] * R[0] * (0.5f * k + R[0] * (parm_a + parm_b * R[0]));
    dV[0]  = R[0] * (k + R[0] * (3 * parm_a + 4 * parm_b * R[0]));
    ddV[0] = k + 6 * R[0] * (parm_a + 2 * parm_b * R[0]);
    return 0;
}

int SCTEST_ForceField::ForceField_npes_SC2D(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                            const int& rdim) {
    const double k0     = mod_M[1] * mod_W[1] * mod_M[1];  // = 2.0f
    const double k1     = mod_M[1] * mod_W[1] * mod_M[1];  // = 1/9
    const double parm_a = -0.1f, parm_b = 0.1f, parm_c = 2.0f;

    V[0]  = R[0] * R[0] * (0.5 * k0 + R[0] * (parm_a + parm_b * R[0])) + 0.5f * k1 * R[1] * R[1] + parm_c * R[0] * R[1];
    dV[0] = R[0] * (k0 + R[0] * (3 * parm_a + 4 * parm_b * R[0])) + parm_c * R[1];
    dV[1] = k1 * R[1] + parm_c * R[0];
    ddV[0] = k0 + 6 * R[0] * (parm_a + 2 * parm_b * R[0]);
    ddV[1] = ddV[2] = parm_c;
    ddV[3]          = k1;
    return 0;
}

int SCTEST_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                       const int& rdim) {
    switch (fftype) {
        case SCTESTPolicy::SC1D: {
            ForceField_npes_SC1D(V, dV, ddV, R, P, flag, rdim);
            break;
        }
        case SCTESTPolicy::SC2D: {
            ForceField_npes_SC2D(V, dV, ddV, R, P, flag, rdim);
            break;
        }
    }
    return 0;
}

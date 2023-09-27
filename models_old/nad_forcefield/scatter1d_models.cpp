#include "scatter1d_models.h"

#include "../../utils/definitions.h"

Scatter1D_ForceField::Scatter1D_ForceField(const Param& iparm, const int& child)
    : Nad_ForceField(iparm){};  // for child to call

Scatter1D_ForceField::Scatter1D_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    CHECK_EQ(N, 1);  // N must be 1
    CHECK_EQ(F, 2);  // F must be 2

    double m0     = Param_GetV(m0, parm, 2000.0f);
    double r0     = Param_GetV(r0, parm, 100.0f);
    double p0     = Param_GetV(p0, parm, 100.0f);
    double varr   = Param_GetV(varr, parm, 0.5f);
    double varp   = Param_GetV(varp, parm, 0.5f);
    mod_M[0]      = m0;
    mod_R0[0]     = r0;
    mod_P0[0]     = p0;
    mod_sigmaR[0] = sqrt(varr);
    mod_sigmaP[0] = sqrt(varp);
    spec          = Param_GetV(spec, parm, 1);

    try {
        Param_GetV(ffflag, parm, "sac");
        fftype = Scatter1DPolicy::_dict.at(ffflag);
    } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }

    switch (fftype) {
        case Scatter1DPolicy::SAC:
        case Scatter1DPolicy::SAC2:
            Xrange = 5.0f;
            break;
        case Scatter1DPolicy::DAC:
            Xrange = 10.0f;
            break;
        case Scatter1DPolicy::ECR:
            Xrange = 13.0f;
            break;
        case Scatter1DPolicy::DBG:
            Xrange = 20.0f;
            break;
        case Scatter1DPolicy::DAG:
            Xrange = 20.0f;
            break;
        case Scatter1DPolicy::DRN:
            Xrange = 10.0f;
            break;
    }

    int plot = Param_GetV(plot, parm, 0);
    if (plot > 0) Scatter1D_plot(Xrange);

    int parm_ctype = Param_GetV(parm_ctype, parm, 0);
    switch (parm_ctype) {
        case 0:  // nothing happens
            break;
        case 1: {
            double parm_i = Param_GetV(parm_i, parm);
            mod_P0[0]     = sqrt(10.0f * parm_i);
            break;
        }
        case -1: {
            double parm_i = Param_GetV(parm_i, parm);
            double Erev   = 10.0f * parm_i / (2 * mod_M[0]) - 0.005f;
            CHECK_GT(Erev, 0);
            mod_P0[0] = sqrt(2 * mod_M[0] * Erev);
            break;
        }
        case -2: {
            double parm_i = Param_GetV(parm_i, parm);
            double Erev   = 10.0f * parm_i / (2 * mod_M[0]) - 0.015f;
            CHECK_GT(Erev, 0);
            mod_P0[0] = sqrt(2 * mod_M[0] * Erev);
            break;
        }
        case 2: {
            double parm_E = Param_GetV(parm_E, parm);
            mod_P0[0]     = sqrt(2 * mod_M[0] * parm_E);
            break;
        }
    }
    int parm_wtype = Param_GetV(parm_wtype, parm, 0);
    switch (parm_wtype) {
        case 0:  // nothing happens
            break;
        case 1: {
            double parm_s = Param_GetV(parm_s, parm);
            mod_sigmaP[0] = sqrt(0.5f) / parm_s;
            break;
        }
        case 2: {
            double parm_a = Param_GetV(parm_a, parm);
            mod_sigmaP[0] = sqrt(0.5f * parm_a);
            break;
        }
    }

    if (std::abs(mod_R0[0]) < Xrange) LOG(WARNING) << "mod_R0 should far away from Xrange";
    if (mod_tend < 0) mod_tend = std::abs(5 * mod_R0[0] * mod_M[0] / mod_P0[0]);
    if (mod_dt < 0) {
        double dt_opt = std::abs(Xrange * mod_M[0] / mod_P0[0]) / 5000;  // 10000 steps
        std::stringstream sstr;
        sstr << FMT(1) << dt_opt;
        sstr >> dt_opt;
        Param_Reset(mod_dt, dt_opt);
    }

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

int Scatter1D_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                          int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    return 0;
}

int Scatter1D_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim, const int& fdim) {
    if (spec == 1) {  // specify forward and backward
        return (rdim > 0 && fdim > 0) ? ((np[0] > 0) ? 0 : 1) : 2;
    }
    return Nad_ForceField ::ForceField_spec(nr, np, nm, rdim, fdim);
}

int Scatter1D_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                          const int& rdim) {
    V[0] = 0.0f;
    if (flag < 1) return 0;
    dV[0] = 0.0f;
    if (flag < 2) return 0;
    ddV[0] = 0.0f;
    return 0;
}

int Scatter1D_ForceField::Scatter1D_plot(const double& Xrange) {
    std::ofstream ofs(std::string("surf_") + ffflag);
    double *R = new double, *V = new double[FF], *dV = new double[NFF], *ddV;
    double *E = new double[F], *T = new double[FF], *dE = new double[NFF], *ddE;
    double* workr      = new double[FF];
    num_complex* workc = new num_complex[FF];

    ofs << FMT(8) << "R";
    for (int i = 0; i < F; ++i)
        for (int j = 0; j < F; ++j) ofs << FMT(8) << utils::concat("V", i, j);
    for (int k = 0; k < N; ++k)
        for (int i = 0; i < F; ++i)
            for (int j = 0; j < F; ++j) ofs << FMT(8) << utils::concat("dV", k, i, j);
    ofs << std::endl;

    for (int i = -100; i <= 100; ++i) {
        R[0] = 0.01 * Xrange * i;
        ForceField_epes(V, dV, ddV, R, 1, N, F);
        ofs << FMT(8) << R[0];
        for (int k = 0; k < FF; ++k) ofs << FMT(8) << V[k];
        for (int k = 0; k < NFF; ++k) ofs << FMT(8) << dV[k];
        ofs << std::endl;
    }
    ofs.close();
    delete R, delete[] V, delete[] dV, delete[] E, delete[] dE, delete[] T, delete[] workr, delete[] workc;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag,
                                          const int& rdim, const int& fdim) {
    switch (fftype) {
        case Scatter1DPolicy::SAC:
            ForceField_epes_SAC(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::SAC2:
            ForceField_epes_SAC2(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::DAC:
            ForceField_epes_DAC(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::ECR:
            ForceField_epes_ECR(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::DBG:
            ForceField_epes_DBG(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::DAG:
            ForceField_epes_DAG(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case Scatter1DPolicy::DRN:
            ForceField_epes_DRN(V, dV, ddV, R, flag, rdim, fdim);
            break;
        default:
            LOG(FATAL);
    }
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_SAC(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double V0 = 0.01f, V1 = 0.005f, E0 = 0.0f, a = 1.6f, b = 1.0f;

    V[0] = (R[0] > 0) ? V0 * (1.0 - exp(-a * R[0])) : -V0 * (1.0 - exp(a * R[0]));
    V[3] = -V[0];
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = (R[0] > 0) ? V0 * a * exp(-a * R[0]) : V0 * a * exp(a * R[0]);
    dV[3] = -dV[0];
    dV[1] = -2 * b * R[0] * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;
    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_SAC2(double* V, double* dV, double* ddV, double* R, const int& flag,
                                               const int& rdim, const int& fdim) {
    const double V0 = 0.01f, V1 = 0.005f, E0 = 0.0f, a = 1.6f, b = 1.0f;

    V[0] = V0 * tanh(a * R[0]);
    V[3] = -V[0];
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    double tmp = cosh(a * R[0]);
    dV[0]      = V0 * a / (tmp * tmp);
    dV[3]      = -dV[0];
    dV[1]      = -2 * b * R[0] * V[1];
    dV[2]      = dV[1];
    if (flag < 2) return 0;
    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}


int Scatter1D_ForceField::ForceField_epes_DAC(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double V0 = 0.10f, V1 = 0.015f, E0 = 0.05f, a = 0.28f, b = 0.06f;

    V[0] = 0.0f;
    V[3] = -V0 * exp(-a * R[0] * R[0]) + E0;
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0f;
    dV[3] = 2 * a * R[0] * V0 * exp(-a * R[0] * R[0]);
    dV[1] = -2 * b * R[0] * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_ECR(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double V0 = 6.0e-4, V1 = 0.1f, E0 = 0.0f, a = 0.0f, b = 0.9f;

    V[0] = -V0, V[3] = V0;
    V[1] = (R[0] < 0) ? V1 * exp(b * R[0]) : V1 * (2.0 - exp(-b * R[0]));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0f, dV[3] = 0.0f;
    dV[1] = (R[0] < 0) ? V1 * b * exp(b * R[0]) : V1 * b * exp(-b * R[0]);
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_DBG(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double V0 = 6.0e-4, V1 = 0.1f, E0 = 0.0f, a = 0.0f, b = 0.9f, Z = 10.0f;

    V[0] = -V0, V[3] = V0;
    if (R[0] < -Z) {
        V[1] = V1 * exp(b * (R[0] - Z)) + V1 * (2.0 - exp(b * (R[0] - Z)));
    } else if (R[0] < Z) {
        V[1] = V1 * exp(b * (R[0] - Z)) + V1 * exp(-b * (R[0] + Z));
    } else {
        V[1] = V1 * exp(-b * (R[0] + Z)) + V1 * (2.0 - exp(-b * (R[0] - Z)));
    }
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0f, dV[3] = 0.0f;
    if (R[0] < -Z) {
        dV[1] = V1 * exp(b * (R[0] - Z)) * b + V1 * (-exp(b * (R[0] - Z))) * b;
    } else if (R[0] < Z) {
        dV[1] = V1 * exp(b * (R[0] - Z)) * b + V1 * (exp(-b * (R[0] + Z))) * (-b);
    } else {
        dV[1] = V1 * exp(-b * (R[0] + Z)) * (-b) + V1 * (-exp(-b * (R[0] - Z))) * (-b);
    }
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_DAG(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double V0 = 6.0e-4, V1 = 0.1f, E0 = 0.0f, a = 0.0f, b = 0.9f, Z = 10.0f;

    V[0] = -V0, V[3] = V0;
    if (R[0] < -Z) {
        V[1] = -V1 * exp(b * (R[0] - Z)) + V1 * exp(b * (R[0] + Z));
    } else if (R[0] < Z) {
        V[1] = -V1 * exp(b * (R[0] - Z)) - V1 * exp(-b * (R[0] + Z)) + 2 * V1;
    } else {
        V[1] = V1 * exp(-b * (R[0] - Z)) - V1 * exp(-b * (R[0] + Z));
    }
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0f, dV[3] = 0.0f;
    if (R[0] < -Z) {
        dV[1] = -V1 * exp(b * (R[0] - Z)) * b + V1 * exp(b * (R[0] + Z)) * b;
    } else if (R[0] < Z) {
        dV[1] = -V1 * exp(b * (R[0] - Z)) * b - V1 * exp(-b * (R[0] + Z)) * (-b);
    } else {
        dV[1] = V1 * exp(-b * (R[0] - Z)) * (-b) - V1 * exp(-b * (R[0] + Z)) * (-b);
    }
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int Scatter1D_ForceField::ForceField_epes_DRN(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    const double E0 = 0.01f, V1 = 0.03f, b = 3.2f, Z = 2.0f;

    V[0] = 0.0f;
    V[3] = E0;
    V[1] = V1 * (exp(-b * (R[0] - Z) * (R[0] - Z)) + exp(-b * (R[0] + Z) * (R[0] + Z)));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0f;
    dV[3] = 0.0f;
    dV[1] = V1 * (exp(-b * (R[0] - Z) * (R[0] - Z)) * (-2 * b * (R[0] - Z)) +
                  exp(-b * (R[0] + Z) * (R[0] + Z)) * (-2 * b * (R[0] + Z)));
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

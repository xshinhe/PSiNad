#include "nad1d_models.h"

#include "../../utils/definitions.h"

NAD1D_ForceField::NAD1D_ForceField(const Param& iparm, const int& child)
    : Nad_ForceField(iparm){};  // for child to call
NAD1D_ForceField::NAD1D_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    CHECK_EQ(N, 1);  // N must be 1

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
    spec          = Param_GetV(spec, parm, 0);

    try {
        Param_GetV(ffflag, parm, "morse3a");
        fftype = NAD1DPolicy::_dict.at(ffflag);
    } catch (std::runtime_error& e) { LOG(FATAL) << e.what(); }

    switch (fftype) {
        case NAD1DPolicy::MORSE3A: {
            CHECK_EQ(F, 3);
            mod_M[0]       = 20000.0f;
            mod_R0[0]      = 2.9f;
            mod_P0[0]      = 0.0f;
            double wground = 5.0e-03;
            mod_sigmaR[0]  = sqrt(0.5f / (mod_M[0] * wground));
            mod_sigmaP[0]  = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::MORSE3B: {
            CHECK_EQ(F, 3);
            mod_M[0]       = 20000.0f;
            mod_R0[0]      = 3.3f;
            mod_P0[0]      = 0.0f;
            double wground = 5.0e-03;
            mod_sigmaR[0]  = sqrt(0.5f / (mod_M[0] * wground));
            mod_sigmaP[0]  = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::MORSE3C: {
            CHECK_EQ(F, 3);
            mod_M[0]       = 20000.0f;
            mod_R0[0]      = 2.1f;
            mod_P0[0]      = 0.0f;
            double wground = 5.0e-03;
            mod_sigmaR[0]  = sqrt(0.5f / (mod_M[0] * wground));
            mod_sigmaP[0]  = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::MORSE15: {
            CHECK_EQ(F, 15);
            mod_M[0]  = 2000.0f;
            mod_R0[0] = 13.0f;
            mod_P0[0] = -30.0f;
            // 2*alpha^2 * D = m*w^2
            const double Dg = 0.2f, alpha = 0.4f;
            mod_sigmaR[0] = sqrt(0.5f / std::sqrt(mod_M[0] * 2 * alpha * alpha * Dg));
            mod_sigmaP[0] = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::IVP1:
        case NAD1DPolicy::IVP2: {
            CHECK_EQ(F, 2);
            double m0     = Param_GetV(m0, parm, 100.0f);
            double w0     = Param_GetV(w0, parm, 0.03f);
            mod_M[0]      = m0;
            mod_W[0]      = w0;
            mod_R0[0]     = 0.0f;
            mod_P0[0]     = 10.0f;
            mod_sigmaR[0] = sqrt(0.5f / (mod_M[0] * mod_W[0]));
            mod_sigmaP[0] = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::IVP3: {
            CHECK_EQ(F, 2);
            double m0     = Param_GetV(m0, parm, 100.0f);
            double w0     = Param_GetV(w0, parm, 0.03f);
            mod_M[0]      = m0;
            mod_W[0]      = w0;
            mod_R0[0]     = -0.8f;
            mod_P0[0]     = 0.0f;
            mod_sigmaR[0] = sqrt(0.5f / (mod_M[0] * mod_W[0]));
            mod_sigmaP[0] = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::IVP4: {
            CHECK_EQ(F, 2);
            double m0     = Param_GetV(m0, parm, 100.0f);
            double w0     = Param_GetV(w0, parm, 0.03f);
            mod_M[0]      = m0;
            mod_W[0]      = w0;
            mod_R0[0]     = -0.8f;
            mod_P0[0]     = 0.0f;
            mod_sigmaR[0] = sqrt(0.5f / (mod_M[0] * mod_W[0]));
            mod_sigmaP[0] = 0.5f / mod_sigmaR[0];
            break;
        }
        case NAD1DPolicy::CL1D: {
            CHECK_EQ(F, 2);
            double cl1d_e = Param_GetV(cl1d_e, parm, 1.0f);
            pm[0]         = cl1d_e;
            double cl1d_d = Param_GetV(cl1d_d, parm, 1.0f);
            pm[1]         = cl1d_d;
            double cl1d_c = Param_GetV(cl1d_c, parm, 1.0f);
            pm[2]         = cl1d_c;
            double cl1d_w = Param_GetV(cl1d_w, parm, 1.0f);
            mod_W[0]      = cl1d_w;
            mod_M[0]      = 1.0f;
            break;
        }
        case NAD1DPolicy::JC1D: {
            CHECK_EQ(F, 2);
            double jc1d_e = Param_GetV(jc1d_e, parm, 1.0f);
            pm[0]         = jc1d_e;
            double jc1d_d = Param_GetV(jc1d_d, parm, 1.0f);
            pm[1]         = jc1d_d;
            double jc1d_c = Param_GetV(jc1d_c, parm, 1.0f);
            pm[2]         = jc1d_c;
            double jc1d_w = Param_GetV(jc1d_w, parm, 1.0f);
            mod_W[0]      = jc1d_w;
            mod_M[0]      = 1.0f;
            break;
        }
        case NAD1DPolicy::NA_I: {
            CHECK_EQ(F, 2);
            double parm_E = Param_GetV(parm_E, parm, 1.0f);
            pm[0]         = parm_E;
            pm[1]         = 0.0f;  // l^2
            // initializetion in au
            mod_M[0]  = (22.989769282f * 126.904473f) / (22.989769282f + 126.904473f) / phys::au_2_amu;
            mod_R0[0] = 160.0f / phys::au_2_ang;
            mod_P0[0] = -sqrt(2 * mod_M[0] * pm[0]);

            double vel = sqrt(2 * parm_E / mod_M[0]);
            if (mod_tend < 0) Param_Reset(mod_tend, 2 * mod_R0[0] / vel);
            if (mod_dt < 0) Param_Reset(mod_dt, mod_tend / 50000);
            break;
        }
    }

    int plot = Param_GetV(plot, parm, 0);
    if (plot > 0) {
        NAD1D_plot();
        exit(0);
    }

    if (mod_tend < 0) Param_Reset(mod_tend, 80.0f);  // for IVP1-4, run 1000 au
    if (mod_dt < 0) Param_Reset(mod_dt, 0.02f);

    tag = name() + ffflag + "_" + tag;
    CheckForceField();
};

int NAD1D_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                      int& eocc, const int& rdim, const int& fdim, const int& icycle) {
    switch (fftype) {
        case NAD1DPolicy::NA_I: {
            const double bmax = 9.0f / phys::au_2_ang;
            double randu;
            rand_uniform(&randu);
            // randu = 0.5f;
            double b = sqrt(randu) * bmax;
            pm[1]    = 2 * mod_M[0] * pm[0] * b * b;  // l^2
            nm[0]    = mod_M[0];
            nr[0]    = mod_R0[0];
            np[0]    = mod_P0[0];
            break;
        }
        default:
            Nad_ForceField::ForceField_init_default(nr, np, nm, erho, eeac, eocc, rdim, fdim, icycle);
    }
    return 0;
}

int NAD1D_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim, const int& fdim) {
    if (spec == 1) {  // specify forward and backward
        return (rdim > 0 && fdim > 0) ? ((np[0] > 0) ? 0 : 1) : 2;
    }
    return Nad_ForceField ::ForceField_spec(nr, np, nm, rdim, fdim);
}

int NAD1D_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                      const int& rdim) {
    V[0] = 0.0f;
    if (flag < 1) return 0;
    dV[0] = 0.0f;
    if (flag < 2) return 0;
    ddV[0] = 0.0f;
    return 0;
}

int NAD1D_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                      const int& fdim) {
    switch (fftype) {
        case NAD1DPolicy::MORSE3A:
            ForceField_epes_Morse3A(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::MORSE3B:
            ForceField_epes_Morse3B(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::MORSE3C:
            ForceField_epes_Morse3C(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::MORSE15:
            ForceField_epes_Morse15(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::IVP1:
            ForceField_epes_IVP1(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::IVP2:
            ForceField_epes_IVP2(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::IVP3:
            ForceField_epes_IVP3(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::IVP4:
            ForceField_epes_IVP4(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::CL1D:
            ForceField_epes_CL1D(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::JC1D:
            ForceField_epes_JC1D(V, dV, ddV, R, flag, rdim, fdim);
            break;
        case NAD1DPolicy::NA_I:
            ForceField_epes_NA_I(V, dV, ddV, R, flag, rdim, fdim);
            break;
    }
    return 0;
}

int NAD1D_ForceField::NAD1D_plot() {
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

    for (int i = 0; i <= 1000; ++i) {
        R[0] = 2 + 0.001 * 20 * i;
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

int NAD1D_ForceField::ForceField_epes_Morse3A(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    // V
    const double De[3]    = {0.003f, 0.004f, 0.003f};
    const double Be[3]    = {0.65f, 0.60f, 0.65f};
    const double Re[3]    = {5.0f, 4.0f, 6.0f};
    const double C[3]     = {0.0f, 0.01f, 0.006f};
    const double Aijx[3]  = {0.002f, 0.00f, 0.002f};
    const double Aeijx[3] = {16.0f, 0.00f, 16.0f};
    const double Rijx[3]  = {4.80f, 0.00f, 3.40f};

    V[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) + C[0];
    V[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) + C[1];
    V[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * exp(-Be[0] * (R[0] - Re[0])) * 2 * Be[0];
    dV[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * exp(-Be[1] * (R[0] - Re[1])) * 2 * Be[1];
    dV[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * exp(-Be[2] * (R[0] - Re[2])) * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_Morse3B(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    // V
    const double De[3]    = {0.020f, 0.010f, 0.003f};
    const double Be[3]    = {0.65f, 0.40f, 0.65f};
    const double Re[3]    = {4.5f, 4.0f, 4.4f};
    const double C[3]     = {0.0f, 0.01f, 0.02f};
    const double Aijx[3]  = {0.00f, 0.005f, 0.005f};
    const double Aeijx[3] = {0.00f, 32.0f, 32.0f};
    const double Rijx[3]  = {0.00f, 3.34f, 3.66f};

    V[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) + C[0];
    V[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) + C[1];
    V[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * exp(-Be[0] * (R[0] - Re[0])) * 2 * Be[0];
    dV[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * exp(-Be[1] * (R[0] - Re[1])) * 2 * Be[1];
    dV[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * exp(-Be[2] * (R[0] - Re[2])) * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_Morse3C(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    // V
    const double De[3]    = {0.020f, 0.020f, 0.003f};
    const double Be[3]    = {0.40f, 0.65f, 0.65f};
    const double Re[3]    = {4.0f, 4.5f, 6.0f};
    const double C[3]     = {0.02f, 0.00f, 0.02f};
    const double Aijx[3]  = {0.00f, 0.005f, 0.005f};
    const double Aeijx[3] = {0.00f, 32.0f, 32.0f};
    const double Rijx[3]  = {0.00f, 4.97f, 3.40f};

    V[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) + C[0];
    V[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) + C[1];
    V[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0f - exp(-Be[0] * (R[0] - Re[0]))) * exp(-Be[0] * (R[0] - Re[0])) * 2 * Be[0];
    dV[4] = De[1] * (1.0f - exp(-Be[1] * (R[0] - Re[1]))) * exp(-Be[1] * (R[0] - Re[1])) * 2 * Be[1];
    dV[8] = De[2] * (1.0f - exp(-Be[2] * (R[0] - Re[2]))) * exp(-Be[2] * (R[0] - Re[2])) * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}


int NAD1D_ForceField::ForceField_epes_Morse15(double* V, double* dV, double* ddV, double* R, const int& flag,
                                              const int& rdim, const int& fdim) {
    // V
    const double Dg = 0.2f, De = 0.05f, Dc = 0.05f, alpha = 0.4f, lambda = 1.0f, eta = 0.004f;
    const double Re[15] = {3.0000000000000f,    6.8335793401926175f, 6.909940519903887f, 6.988372712202933f,
                           7.0690037103879675f, 7.151973244334301f,  7.237434522736795f, 7.325556032471846f,
                           7.416523648311893f,  7.5105431195440095f, 7.607843017311827f, 7.708678249086644f,
                           7.813334276495011f,  7.922132212503321f,  8.035435027584366};

    // {3.0000000000000000f, 0.9146012197794767f, 0.9044940375386811f, 0.8943426828231982f, 0.8841415645058391f,
    //  0.8738847009741216f, 0.8635656710074937f, 0.853177557074951f,  0.8427128795608424f, 0.8321635200703644f,
    //  0.8215206315085599f, 0.8107745320334944f, 0.7999145792087743f, 0.7889290196565726f, 0.7778048081460116f};

    for (int i = 0; i < FF; ++i) V[i] = 0.0f;
    V[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * (1 - exp(-alpha * (R[0] - Re[0])));
    for (int i = 1; i < F; ++i) {
        V[i * (F + 1)] = exp(-alpha * R[0]) + eta * (i + 1) + De;
        V[i]           = Dc * exp(-lambda * (R[0] - Re[i]) * (R[0] - Re[i]));
        V[i * F]       = V[i];
    }
    if (flag < 1) return 0;

    for (int i = 0; i < NFF; ++i) dV[i] = 0.0f;
    dV[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * exp(-alpha * (R[0] - Re[0])) * 2 * alpha;
    for (int i = 1; i < F; ++i) {
        dV[i * (F + 1)] = -alpha * exp(-alpha * R[0]);
        dV[i]           = V[i] * (-2 * lambda * (R[0] - Re[i]));
        dV[i * F]       = dV[i];
    }

    if (flag < 2) return 0;
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_IVP1(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double A = 0.003f, B = 0.004f, C = 0.05f,  // @DIFF
        a = 0.4f, b = 1.0f, Re = 4.0f, Rc = 1.5f;

    // @deps. on mod_M and mod_W
    V[0] = 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * R[0] * R[0];
    V[3] = A * (1.0f - exp(-a * (R[0] - Re))) * (1.0f - exp(-a * (R[0] - Re))) + B;
    V[1] = C * (1.0f - std::tanh(b * (R[0] - Rc)));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = mod_M[0] * mod_W[0] * mod_W[0] * R[0];
    dV[3] = A * (1.0f - exp(-a * (R[0] - Re))) * exp(-a * (R[0] - Re)) * 2 * a;
    dV[1] = C * b / (std::cosh(b * (R[0] - Rc)) * std::cosh(b * (R[0] - Rc)));
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_IVP2(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double A = 0.003f, B = 0.004f, C = 0.01f,  // @DIFF
        a = 0.4f, b = 1.0f, Re = 4.0f, Rc = 1.5f;

    // @deps. on mod_M and mod_W
    V[0] = 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * R[0] * R[0];
    V[3] = A * (1.0f - exp(-a * (R[0] - Re))) * (1.0f - exp(-a * (R[0] - Re))) + B;
    V[1] = C * (1.0f - std::tanh(b * (R[0] - Rc)));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = mod_M[0] * mod_W[0] * mod_W[0] * R[0];
    dV[3] = A * (1.0f - exp(-a * (R[0] - Re))) * exp(-a * (R[0] - Re)) * 2 * a;
    dV[1] = C * b / (std::cosh(b * (R[0] - Rc)) * std::cosh(b * (R[0] - Rc)));
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_IVP3(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double A = 0.0001f, B = 0.1f, C = 0.1f,  // @DIFF
        Re             = 0.8f;
    const double alpha = 0.5f;
    const double a     = 0.8f;  // mod_W[0] * sqrt(mod_M[0] / (2 * A));

    // @deps. on mod_M and mod_W
    V[0] = alpha * 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * (R[0] + Re) * (R[0] + Re) +
           (1 - alpha) * A * (1.0f - exp(-a * (R[0] + Re))) * (1.0f - exp(-a * (R[0] + Re))) + B;
    V[3] = alpha * (0.5f * mod_M[0] * mod_W[0] * mod_W[0] * (R[0] - Re) * (R[0] - Re)) +
           (1 - alpha) * A * (1.0f - exp(-a * (R[0] - Re))) * (1.0f - exp(-a * (R[0] - Re))) - B;
    V[1] = C;
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = alpha * (mod_M[0] * mod_W[0] * mod_W[0] * (R[0] + Re)) +
            (1 - alpha) * A * (1.0f - exp(-a * (R[0] + Re))) * exp(-a * (R[0] + Re)) * 2 * a;
    dV[3] = alpha * (mod_M[0] * mod_W[0] * mod_W[0] * (R[0] - Re)) +
            (1 - alpha) * A * (1.0f - exp(-a * (R[0] - Re))) * exp(-a * (R[0] - Re)) * 2 * a;
    dV[1] = 0;
    dV[2] = 0;
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_IVP4(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double A = 0.00001f, B = 0.005f, C = 0.01f,  // @DIFF
        Re             = 0.8f;
    const double alpha = 0.5f;
    const double a     = 0.8f;  // mod_W[0] * sqrt(mod_M[0] / (2 * A));

    // @deps. on mod_M and mod_W
    V[0] = alpha * 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * (R[0] + Re) * (R[0] + Re) +
           (1 - alpha) * A * (1.0f - exp(-a * (R[0] + Re))) * (1.0f - exp(-a * (R[0] + Re))) + B;
    V[3] = alpha * (0.5f * mod_M[0] * mod_W[0] * mod_W[0] * (R[0] - Re) * (R[0] - Re)) +
           (1 - alpha) * A * (1.0f - exp(-a * (R[0] - Re))) * (1.0f - exp(-a * (R[0] - Re))) - B;
    V[1] = C;
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = alpha * (mod_M[0] * mod_W[0] * mod_W[0] * (R[0] + Re)) +
            (1 - alpha) * A * (1.0f - exp(-a * (R[0] + Re))) * exp(-a * (R[0] + Re)) * 2 * a;
    dV[3] = alpha * (mod_M[0] * mod_W[0] * mod_W[0] * (R[0] - Re)) +
            (1 - alpha) * A * (1.0f - exp(-a * (R[0] - Re))) * exp(-a * (R[0] - Re)) * 2 * a;
    dV[1] = 0;
    dV[2] = 0;
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_CL1D(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    // @deps. on mod_M and mod_W
    double meanV = 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * R[0] * R[0];
    V[0]         = meanV + pm[0] + pm[2] * R[0];
    V[3]         = meanV - pm[0] - pm[2] * R[0];
    V[1]         = pm[1];
    V[2]         = pm[1];
    if (flag < 1) return 0;

    dV[0] = mod_M[0] * mod_W[0] * mod_W[0] * R[0] + pm[2];
    dV[3] = -dV[0];
    dV[1] = 0;
    dV[2] = 0;
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_JC1D(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    // @deps. on mod_M and mod_W
    double meanV = 0.5f * mod_M[0] * mod_W[0] * mod_W[0] * R[0] * R[0];
    V[0]         = meanV + pm[0];
    V[3]         = meanV - pm[0];
    V[1]         = pm[1] + pm[2] * R[0];
    V[2]         = pm[1] + pm[2] * R[0];
    if (flag < 1) return 0;

    dV[0] = mod_M[0] * mod_W[0] * mod_W[0] * R[0];
    dV[3] = -dV[0];
    dV[1] = pm[2];
    dV[2] = pm[2];
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

int NAD1D_ForceField::ForceField_epes_NA_I(double* V, double* dV, double* ddV, double* R, const int& flag,
                                           const int& rdim, const int& fdim) {
    const double au_2_ev  = 27.21138602f;
    const double au_2_ang = 0.529177249f;
    const double Acov     = 3150.0f / au_2_ev;
    const double Bcov12   = pow(2.647f / au_2_ang, 12) / au_2_ev;
    const double Ccov     = 1000.0f / au_2_ev / pow(au_2_ang, 6);
    const double Rcov     = 0.435f / au_2_ang;
    const double Aion     = 2760.0f / au_2_ev;
    const double Bion8    = pow(2.398f / au_2_ang, 8) / au_2_ev;
    const double Cion     = 11.3f / au_2_ev / pow(au_2_ang, 6);
    const double Rion     = 0.3489f / au_2_ang;
    const double al_Mp    = 0.408f / pow(au_2_ang, 3);
    const double al_Xm    = 6.431f / pow(au_2_ang, 3);
    const double al_M     = 27.0f / pow(au_2_ang, 3);
    const double al_X     = 7.0f / pow(au_2_ang, 3);
    const double Eth      = 2.075f / au_2_ev;
    const double Aoff     = 17.08f / au_2_ev;
    const double Roff     = 1.239f / au_2_ang;  // @bugs?

    double r   = R[0];
    double r2  = r * r;
    double r4  = r2 * r2;
    double r5  = r4 * r;
    double r6  = r2 * r4;
    double r7  = r6 * r;
    double r8  = r4 * r4;
    double r9  = r8 * r;
    double r12 = r6 * r6;
    double r13 = r12 * r;

    double Vcent = pm[1] / (2 * mod_M[0] * r2);

    V[0] = Vcent                                     //
           + (Acov + Bcov12 / r12) * exp(-r / Rcov)  //
           - Ccov / r6;                              //
    V[3] = Vcent                                     //
           + (Aion + Bion8 / r8) * exp(-r / Rion)    //
           - 1 / r                                   //
           - 0.5 * (al_Mp + al_Xm) / r4              //
           - Cion / r6                               //
           - 2 * al_Mp * al_Xm / r7                  //
           + Eth;
    V[1] = Aoff * exp(-r / Roff);
    V[2] = V[1];

    if (flag < 1) return 0;

    double dVcent = -2 * pm[1] / (2 * mod_M[0] * r2 * r);

    dV[0] = dVcent                                                  //
            + (-12 * Bcov12 / r13) * exp(-r / Rcov)                 //
            + (Acov + Bcov12 / r12) * (-1 / Rcov) * exp(-r / Rcov)  //
            + 6 * Ccov / r7;                                        //
    dV[3] = dVcent                                                  //
            + (-8 * Bion8 / r9) * exp(-r / Rion)                    //
            + (Aion + Bion8 / r8) * (-1 / Rion) * exp(-r / Rion)    //
            + 1 / r2                                                //
            + 2 * (al_Mp + al_Xm) / r5                              //
            + 6 * Cion / r7                                         //
            + 14 * al_Mp * al_Xm / r8;                              //
    dV[1] = (-1 / Roff) * Aoff * exp(-r / Roff);
    dV[2] = dV[1];

    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < NNFF; ++i) ddV[i] = 0;
    return 0;
}

#include "Model_NAD1D.h"

#include "../kernels/Kernel_Declare.h"
#include "../kernels/Kernel_Random.h"

namespace PROJECT_NS {

double mspes_parm[100];

int mspes_SAC(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_SAC2(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_SAC3(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V1 = 0.04f, V2 = 0.01f, Vc = 0.005f, a = 1.0f, b = 1.0f, Rc = 0.7;

    V[0] = V1 * (1.0e0 + tanh(a * R[0]));
    V[3] = V2 * (1.0e0 - tanh(a * R[0]));
    V[1] = Vc * exp(-b * (R[0] + Rc) * (R[0] + Rc));
    V[2] = V[1];
    if (flag < 1) return 0;

    double tmp = cosh(a * R[0]);
    dV[0]      = +V1 * a / (tmp * tmp);
    dV[3]      = -V2 * a / (tmp * tmp);
    dV[1]      = -2 * b * (R[0] + Rc) * V[1];
    dV[2]      = dV[1];
    if (flag < 2) return 0;

    // ddV
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}


int mspes_DAC(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_ECR(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_DBG(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_DAG(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_DRN(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}


int mspes_MORSE3A(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 9; ++i) ddV[i] = 0;
    return 0;
}

int mspes_MORSE3B(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 9; ++i) ddV[i] = 0;
    return 0;
}

int mspes_MORSE3C(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    // for (int i = 0; i < 9; ++i) ddV[i] = 0;
    return 0;
}


int mspes_MORSE15(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double Dg = 0.2f, De = 0.05f, Dc = 0.05f, alpha = 0.4f, lambda = 1.0f, eta = 0.004f;
    const double Re[15] = {3.0000000000000f,    6.8335793401926175f, 6.909940519903887f, 6.988372712202933f,
                           7.0690037103879675f, 7.151973244334301f,  7.237434522736795f, 7.325556032471846f,
                           7.416523648311893f,  7.5105431195440095f, 7.607843017311827f, 7.708678249086644f,
                           7.813334276495011f,  7.922132212503321f,  8.035435027584366};
    for (int i = 0; i < 225; ++i) V[i] = 0.0f;
    V[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * (1 - exp(-alpha * (R[0] - Re[0])));
    for (int i = 1; i < 15; ++i) {
        V[i * 16] = exp(-alpha * R[0]) + eta * (i + 1) + De;
        V[i]      = Dc * exp(-lambda * (R[0] - Re[i]) * (R[0] - Re[i]));
        V[i * 15] = V[i];
    }
    if (flag < 1) return 0;

    for (int i = 0; i < 225; ++i) dV[i] = 0.0f;
    dV[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * exp(-alpha * (R[0] - Re[0])) * 2 * alpha;
    for (int i = 1; i < 15; ++i) {
        dV[i * 16] = -alpha * exp(-alpha * R[0]);
        dV[i]      = V[i] * (-2 * lambda * (R[0] - Re[i]));
        dV[i * 15] = dV[i];
    }

    if (flag < 2) return 0;
    // for (int i = 0; i < 225; ++i) ddV[i] = 0;
    return 0;
}

int mspes_CL1D(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // @deps. on mass and mod_W
    double meanV = 0.5f * mspes_parm[3] * mspes_parm[3] * R[0] * R[0];  // mass = 1
    V[0]         = meanV + mspes_parm[0] + mspes_parm[2] * R[0];
    V[3]         = meanV - mspes_parm[0] - mspes_parm[2] * R[0];
    V[1]         = mspes_parm[1];
    V[2]         = mspes_parm[1];
    if (flag < 1) return 0;

    dV[0] = mspes_parm[3] * mspes_parm[3] * R[0] + mspes_parm[2];
    dV[3] = -dV[0];
    dV[1] = 0;
    dV[2] = 0;
    if (flag < 2) return 0;

    // ddV
    for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_JC1D(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // @deps. on mass and mod_W
    double meanV = 0.5f * mspes_parm[3] * mspes_parm[3] * R[0] * R[0];  // mass = 1
    V[0]         = meanV + mspes_parm[0];
    V[3]         = meanV - mspes_parm[0];
    V[1]         = mspes_parm[1] + mspes_parm[2] * R[0];
    V[2]         = mspes_parm[1] + mspes_parm[2] * R[0];
    if (flag < 1) return 0;

    dV[0] = mspes_parm[3] * mspes_parm[3] * R[0];
    dV[3] = -dV[0];
    dV[1] = mspes_parm[2];
    dV[2] = mspes_parm[2];
    if (flag < 2) return 0;

    // ddV
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

int mspes_NA_I(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
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
    const double Roff     = 1.239f / au_2_ang;
    const double mred     = (22.989769282f * 126.904473f) / (22.989769282f + 126.904473f) / phys::au_2_amu;

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

    double Vcent = mspes_parm[1] / (2 * mred * r2);

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

    double dVcent = -2 * mspes_parm[1] / (2 * mred * r2 * r);

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
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}


void Model_NAD1D::read_param_impl(Param* PM) {
    nad1d_type = NAD1DPolicy::_from(_Param->get<std::string>("nad1d_flag", LOC(), "SAC"));

    // revise tend and dt
    switch (nad1d_type) {
        case NAD1DPolicy::SAC:
        case NAD1DPolicy::SAC2:
        case NAD1DPolicy::SAC3:
        case NAD1DPolicy::DAC:
        case NAD1DPolicy::ECR:
        case NAD1DPolicy::DBG:
        case NAD1DPolicy::DAG:
        case NAD1DPolicy::DRN: {
            double x0_read   = _Param->get<double>("x0", LOC(), -10.0f);
            double p0_read   = _Param->get<double>("p0", LOC(), 10.0f);
            double mass_read = _Param->get<double>("m0", LOC(), 2000.0f);
            double tend_read = _Param->get<double>("tend", LOC(), -1);
            double dt_read   = _Param->get<double>("dt", LOC(), -1);

            double tend_rev = std::abs(5 * x0_read * mass_read / p0_read);
            double dt_rev   = std::abs(x0_read * mass_read / p0_read) / 5000;
            if (tend_read < 0) (*(_Param->pjson()))["tend"] = tend_rev;
            if (dt_read < 0) (*(_Param->pjson()))["dt"] = dt_rev;
            break;
        }
        case NAD1DPolicy::NA_I: {
            // CHECK_EQ(F, 2);
            // initializetion in au
            double mass_Na = 22.989769282f / phys::au_2_amu;
            double mass_I  = 126.904473f / phys::au_2_amu;
            double mred    = (mass_Na * mass_I) / (mass_Na + mass_I);

            double parm_E    = _Param->get<double>("parm_E", LOC(), 1.0f);
            double tend_read = _Param->get<double>("tend", LOC(), -1);
            double dt_read   = _Param->get<double>("dt", LOC(), -1);

            double x0_max   = 160.0f / phys::au_2_ang;
            double vel      = sqrt(2 * parm_E / mred);
            double tend_rev = 2 * x0_max / vel;
            double dt_rev   = tend_rev / 50000;

            if (tend_read < 0) (*(_Param->pjson()))["tend"] = tend_rev;
            if (dt_read < 0) (*(_Param->pjson()))["dt"] = dt_rev;
            break;
        }
    }
};

void Model_NAD1D::init_data_impl(DataSet* DS) {
    Hsys = DS->reg<num_real>("model.Hsys", Dimension::FF);
    memset(Hsys, 0, Dimension::FF * sizeof(num_real));

    // model field
    mass = DS->reg<double>("model.mass", Dimension::N);
    vpes = DS->reg<double>("model.vpes");                 // not used
    grad = DS->reg<double>("model.grad", Dimension::N);   // not used
    hess = DS->reg<double>("model.hess", Dimension::NN);  // not used
    V    = DS->reg<double>("model.V", Dimension::FF);
    dV   = DS->reg<double>("model.dV", Dimension::NFF);
    // ddV  = DS->reg<double>("model.ddV", Dimension::Dimension::NNFF);

    x0      = DS->reg<double>("model.x0", Dimension::N);
    p0      = DS->reg<double>("model.p0", Dimension::N);
    x_sigma = DS->reg<double>("model.x_sigma", Dimension::N);
    p_sigma = DS->reg<double>("model.p_sigma", Dimension::N);

    // init & integrator
    x      = DS->reg<double>("integrator.x", Dimension::N);
    p      = DS->reg<double>("integrator.p", Dimension::N);
    p_sign = DS->reg<num_complex>("integrator.p_sign", 2);

    double x0_read = _Param->get<double>("x0grid", LOC(), -10.0f);
    int Nxgird     = _Param->get<int>("Nxgrid", LOC(), 101);
    double dx      = (2 * abs(x0_read)) / (Nxgird - 1);
    double* xgrid  = DS->reg<double>("integrator.xgrid", Nxgird);
    for (int i = 0; i < Nxgird; ++i) xgrid[i] = -abs(x0_read) + i * dx;

    DS->reg<double>("init.x", Dimension::N);
    DS->reg<double>("init.p", Dimension::N);

    mass[0]     = _Param->get<double>("m0", LOC(), 2000.0f);
    x0[0]       = _Param->get<double>("x0", LOC(), 100.0f);
    p0[0]       = _Param->get<double>("p0", LOC(), 100.0f);
    double varx = _Param->get<double>("varx", LOC(), 0.5f);
    double varp = _Param->get<double>("varp", LOC(), 0.5f);
    x_sigma[0]  = sqrt(varx);
    p_sigma[0]  = sqrt(varp);

    switch (nad1d_type) {
        case NAD1DPolicy::SAC:
        case NAD1DPolicy::SAC2:
        case NAD1DPolicy::SAC3:  // asymmetrical SAC
        case NAD1DPolicy::DAC:
        case NAD1DPolicy::ECR:
        case NAD1DPolicy::DBG:
        case NAD1DPolicy::DAG:
        case NAD1DPolicy::DRN: {
            mass[0] = 2000.0f;
            break;
        }
        case NAD1DPolicy::MORSE3A:
        case NAD1DPolicy::MORSE3B:
        case NAD1DPolicy::MORSE3C: {
            // CHECK_EQ(F, 3);
            mass[0] = 20000.0f;
            switch (nad1d_type) {
                case NAD1DPolicy::MORSE3A:
                    x0[0] = 2.9f;
                    break;
                case NAD1DPolicy::MORSE3B:
                    x0[0] = 3.3f;
                    break;
                case NAD1DPolicy::MORSE3C:
                    x0[0] = 2.1f;
                    break;
            }
            p0[0]          = 0.0f;
            double wground = 5.0e-03;
            x_sigma[0]     = sqrt(0.5f / (mass[0] * wground));
            p_sigma[0]     = 0.5f / x_sigma[0];
            break;
        }
        case NAD1DPolicy::MORSE15: {
            // CHECK_EQ(F, 15);
            mass[0]   = 2000.0f;
            x0[0]     = 13.0f;
            p0[0]     = -30.0f;
            double Dg = 0.2f, alpha = 0.4f;  // 2*alpha^2 * D = m*w^2
            x_sigma[0] = sqrt(0.5f / std::sqrt(mass[0] * 2 * alpha * alpha * Dg));
            p_sigma[0] = 0.5f / x_sigma[0];
            break;
        }
        case NAD1DPolicy::CL1D:
        case NAD1DPolicy::JC1D: {
            // CHECK_EQ(F, 2);
            mass[0]       = 1.0f;
            mspes_parm[0] = _Param->get<double>("nad1d_e", LOC(), 1.0f);
            mspes_parm[1] = _Param->get<double>("nad1d_d", LOC(), 1.0f);
            mspes_parm[2] = _Param->get<double>("nad1d_c", LOC(), 1.0f);
            mspes_parm[3] = _Param->get<double>("nad1d_w", LOC(), 1.0f);
            break;
        }
        case NAD1DPolicy::NA_I: {
            // CHECK_EQ(F, 2);
            // initializetion in au
            double mass_Na = 22.989769282f / phys::au_2_amu;
            double mass_I  = 126.904473f / phys::au_2_amu;
            mass[0]        = (mass_Na * mass_I) / (mass_Na + mass_I);

            mspes_parm[0] = _Param->get<double>("parm_E", LOC(), 1.0f);
            x0[0]         = 160.0f / phys::au_2_ang;
            p0[0]         = -sqrt(2 * mass[0] * mspes_parm[0]);
            break;
        }
    }
};

void Model_NAD1D::init_calc_impl(int stat) {
    switch (nad1d_type) {
        case NAD1DPolicy::NA_I: {
            double bmax = 9.0f / phys::au_2_ang;
            double randu;
            Kernel_Random::rand_uniform(&randu);
            double b      = sqrt(randu) * bmax;
            mspes_parm[1] = 2 * mass[0] * mspes_parm[0] * b * b;  // l^2 = 2 m E b^2
            for (int j = 0; j < Dimension::N; ++j) {
                x[j] = x0[j];
                p[j] = p0[j];
            }
            break;
        }
        default: {
            Kernel_Random::rand_gaussian(x, Dimension::N);
            Kernel_Random::rand_gaussian(p, Dimension::N);
            for (int j = 0; j < Dimension::N; ++j) {
                x[j] = x0[j] + x[j] * x_sigma[j];
                p[j] = p0[j] + p[j] * p_sigma[j];
            }
            break;
        }
    }
    if (p[0] >= 0) {
        p_sign[0] = phys::math::iu, p_sign[1] = phys::math::iz;
    } else {
        p_sign[0] = phys::math::iz, p_sign[1] = phys::math::iu;
    }
    _DataSet->set("init.x", x, Dimension::N);
    _DataSet->set("init.p", p, Dimension::N);

    exec_kernel(stat);
}

int Model_NAD1D::exec_kernel_impl(int stat) {
    switch (nad1d_type) {
        case NAD1DPolicy::SAC:
            mspes_SAC(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::SAC2:
            mspes_SAC2(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::SAC3:
            mspes_SAC2(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::DAC:
            mspes_DAC(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::ECR:
            mspes_ECR(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::DBG:
            mspes_DBG(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::DAG:
            mspes_DAG(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::DRN:
            mspes_DRN(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3A:
            mspes_MORSE3A(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3B:
            mspes_MORSE3B(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3C:
            mspes_MORSE3C(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::MORSE15:
            mspes_MORSE15(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::JC1D:
            mspes_JC1D(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::CL1D:
            mspes_CL1D(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
        case NAD1DPolicy::NA_I:
            mspes_NA_I(V, dV, ddV, x, 1, 1, Dimension::F);
            break;
    }
    if (p[0] >= 0) {
        p_sign[0] = phys::math::iu, p_sign[1] = phys::math::iz;
    } else {
        p_sign[0] = phys::math::iz, p_sign[1] = phys::math::iu;
    }
    return stat;
}


};  // namespace PROJECT_NS
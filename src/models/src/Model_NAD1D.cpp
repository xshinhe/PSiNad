#include "kids/Model_NAD1D.h"

#include "kids/Kernel_Random.h"
#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"

namespace PROJECT_NS {

double mspes_parm_int[10];
double mspes_parm_real[100];

/**
 * @brief      SAC (simple avoided crossing) model, doi:10.1063/1.459170
 */
int mspes_SAC(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 0.01e0, V1 = 0.005e0, E0 = 0.0e0, a = 1.6e0, b = 1.0e0;

    double exp_aR = (R[0] > 0) ? exp(-a * R[0]) : exp(a * R[0]);
    double dsignR = (R[0] > 0) ? 1.0e0 : -1.0e0;

    V[0] = dsignR * V0 * (1.0e0 - exp_aR);
    V[3] = -V[0];
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = V0 * a * exp_aR;
    dV[3] = -dV[0];
    dV[1] = -2 * b * R[0] * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = -dsignR * V0 * a * a * exp_aR;
    ddV[3] = -dV[0];
    ddV[1] = (-2 * b * V[1] - 2 * b * R[0] * dV[1]);
    ddV[2] = ddV[1];

    return 0;
}

/**
 * @brief      modified SAC (simple avoided crossing) model, doi:10.1063/1.2759932
 */
int mspes_SAC2(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 0.01e0, V1 = 0.005e0, E0 = 0.0e0, a = 1.6e0, b = 1.0e0;

    double tanhaR = tanh(a * R[0]);
    double sechaR = 1.0e0 / cosh(a * R[0]);

    V[0] = V0 * tanhaR;
    V[3] = -V[0];
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = V0 * a * sechaR * sechaR;
    dV[3] = -dV[0];
    dV[1] = -2 * b * R[0] * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = -2 * V0 * a * a * sechaR * sechaR * tanhaR;
    ddV[3] = -dV[0];
    ddV[1] = (-2 * b * V[1] - 2 * b * R[0] * dV[1]);
    ddV[2] = ddV[1];
    return 0;
}

/**
 * @brief      asymmetrical SAC (simple avoided crossing) model, doi:10.1063/1.2759932
 */
int mspes_SAC3(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V1 = 0.04e0, V2 = 0.01e0, Vc = 0.005e0, a = 1.0e0, b = 1.0e0, Rc = 0.7;

    double tanhaR = tanh(a * R[0]);
    double sechaR = 1.0e0 / cosh(a * R[0]);

    V[0] = V1 * (1.0e0 + tanhaR);
    V[3] = V2 * (1.0e0 - tanhaR);
    V[1] = Vc * exp(-b * (R[0] + Rc) * (R[0] + Rc));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = +V1 * a * sechaR * sechaR;
    dV[3] = -V2 * a * sechaR * sechaR;
    dV[1] = -2 * b * (R[0] + Rc) * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = -2 * V1 * a * a * sechaR * sechaR * tanhaR;
    ddV[3] = +2 * V2 * a * a * sechaR * sechaR * tanhaR;
    ddV[1] = (-2 * b * V[1] - 2 * b * (R[0] + Rc) * dV[1]);
    ddV[2] = ddV[1];
    // for (int i = 0; i < 4; ++i) ddV[i] = 0;
    return 0;
}

/**
 * @brief      DAC (dual avoided crossing) model, doi:10.1063/1.459170
 */
int mspes_DAC(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 0.10e0, V1 = 0.015e0, E0 = 0.05e0, a = 0.28e0, b = 0.06e0;

    V[0] = 0.0e0;
    V[3] = -V0 * exp(-a * R[0] * R[0]) + E0;
    V[1] = V1 * exp(-b * R[0] * R[0]);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0e0;
    dV[3] = -2 * a * R[0] * (V[3] - E0);
    dV[1] = -2 * b * R[0] * V[1];
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = 0.0e0;
    ddV[3] = (-2 * a * (V[3] - E0) - 2 * a * R[0] * dV[3]);
    ddV[1] = (-2 * b * V[1] - 2 * b * R[0] * dV[1]);
    ddV[2] = ddV[1];
    return 0;
}

/**
 * @brief      ECR (extended coupling with reflection) model, doi:10.1063/1.459170
 */
int mspes_ECR(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 6.0e-4, V1 = 0.1e0, E0 = 0.0e0, a = 0.0e0, b = 0.9e0;

    double exp_bR = (R[0] > 0) ? exp(-b * R[0]) : exp(b * R[0]);
    double dsignR = (R[0] > 0) ? 1.0e0 : -1.0e0;

    V[0] = -V0, V[3] = V0;
    V[1] = V1 * (1.0e0 + dsignR * (1.0e0 - exp_bR));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0e0, dV[3] = 0.0e0;
    dV[1] = V1 * b * exp_bR;
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    dV[0] = 0.0e0, dV[3] = 0.0e0;
    dV[1] = V1 * (-b * dsignR) * b * exp_bR;
    dV[2] = dV[1];
    return 0;
}

/**
 * @brief      DBG (dumbbell geometry) model,
 * doi:10.1063/1.3506779
 * doi:10.1021/acs.jpclett.0c02533
 */
int mspes_DBG(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 6.0e-4, V1 = 0.1e0, b = 0.9e0, Z = 10.0e0;

    double exp_bR1 = (R[0] < -Z) ? exp(b * (R[0] + Z)) : exp(-b * (R[0] + Z));
    double exp_bR2 = (R[0] < Z) ? exp(b * (R[0] - Z)) : exp(-b * (R[0] - Z));
    double dsignR1 = (R[0] < -Z) ? -1.0e0 : 1.0e0;
    double dsignR2 = (R[0] < Z) ? -1.0e0 : 1.0e0;

    // Piecewise function
    // if (R[0] < -Z) {
    //     V[1] = V1 * (2.0 - exp(b * (R[0] + Z))) + V1 * exp(b * (R[0] - Z));
    // } else if (R[0] < Z) {
    //     V[1] = V1 * exp(-b * (R[0] + Z)) + V1 * exp(b * (R[0] - Z));
    // } else {
    //     V[1] = V1 * exp(-b * (R[0] + Z)) + V1 * (2.0 - exp(-b * (R[0] - Z)));
    // }

    V[0] = -V0, V[3] = V0;
    V[1] = V1 * (1.0e0 - dsignR1 * (1.0e0 - exp_bR1) + 1.0e0 + dsignR2 * (1.0e0 - exp_bR2));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0e0, dV[3] = 0.0e0;
    dV[1] = V1 * b * (-exp_bR1 + exp_bR2);
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = 0.0e0, ddV[3] = 0.0e0;
    ddV[1] = V1 * b * b * (dsignR1 * exp_bR1 - dsignR2 * exp_bR2);
    ddV[2] = ddV[1];
    return 0;
}

/**
 * @brief      DAG (double arch geometry) model,
 * doi:10.1063/1.3506779
 * doi:10.1021/acs.jpclett.0c02533
 */
int mspes_DAG(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double V0 = 6.0e-4, V1 = 0.1e0, b = 0.9e0, Z = 4.0e0;

    double exp_bR1 = (R[0] < -Z) ? exp(b * (R[0] + Z)) : exp(-b * (R[0] + Z));
    double exp_bR2 = (R[0] < Z) ? exp(b * (R[0] - Z)) : exp(-b * (R[0] - Z));
    double dsignR1 = (R[0] < -Z) ? -1.0e0 : 1.0e0;
    double dsignR2 = (R[0] < Z) ? -1.0e0 : 1.0e0;

    // Piecewise function
    // if (R[0] < -Z) {
    //     V[1] = V1 * exp(b * (R[0] + Z)) - V1 * exp(b * (R[0] - Z));
    // } else if (R[0] < Z) {
    //     V[1] = - V1 * exp(-b * (R[0] + Z)) -V1 * exp(b * (R[0] - Z)) + 2 * V1;
    // } else {
    //     V[1] = - V1 * exp(-b * (R[0] + Z)) + V1 * exp(-b * (R[0] - Z)) ;
    // }

    V[0] = -V0, V[3] = V0;
    V[1] = V1 * (1.0e0 + dsignR1 * (1.0e0 - exp_bR1) - 1.0e0 - dsignR2 * (1.0e0 - exp_bR2));
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0e0, dV[3] = 0.0e0;
    dV[1] = V1 * b * (exp_bR1 - exp_bR2);
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = 0.0e0, ddV[3] = 0.0e0;
    ddV[1] = V1 * b * b * (-dsignR1 * exp_bR1 + dsignR2 * exp_bR2);
    ddV[2] = ddV[1];
    return 0;
}

/**
 * @brief      DRN (Rosen-Zener-Demkov noncrossing) model,
 * doi:10.1038/srep24198
 * doi:10.1021/acs.jpclett.0c02533
 */
int mspes_DRN(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double E0 = 0.01e0, V1 = 0.03e0, b = 3.2e0, Z = 2.0e0;

    double exp_bRR1 = exp(-b * (R[0] + Z) * (R[0] + Z));
    double exp_bRR2 = exp(-b * (R[0] - Z) * (R[0] - Z));

    V[0] = 0.0e0;
    V[3] = E0;
    V[1] = V1 * (exp_bRR1 + exp_bRR2);
    V[2] = V[1];
    if (flag < 1) return 0;

    dV[0] = 0.0e0;
    dV[3] = 0.0e0;
    dV[1] = V1 * (exp_bRR1 * (-2 * b * (R[0] + Z)) + exp_bRR2 * (-2 * b * (R[0] - Z)));
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = 0.0e0;
    ddV[3] = 0.0e0;
    ddV[1] = 2 * b * V1 *
             (exp_bRR1 * (2 * b * (R[0] + Z) * (R[0] + Z) - 1)  //
              + exp_bRR2 * (2 * b * (R[0] - Z) * (R[0] - Z) - 1));
    ddV[2] = ddV[1];
    return 0;
}


/**
 * @brief      DPES model,
 * doi:10.1021/acs.jpclett.0c02533
 */
int mspes_DPES(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double E1      = 0.00e0;
    const double E2      = 0.01e0;
    const double E3      = 0.02e0;
    const double A1      = 0.01e0;
    const double A2      = 0.02e0;
    const double B1      = 1.6e0;
    const double B2      = 0.28e0;
    const double A3L[4]  = {0.0025, 0.005, 0.02, 0.04};
    const double B3L[2]  = {0.4, 1.5};
    int          choose1 = mspes_parm_int[0];
    int          choose2 = mspes_parm_int[1];
    int          choose3 = mspes_parm_int[2];
    int          choose4 = mspes_parm_int[3];
    assert(choose1 < 5 && choose2 < 5 && choose3 < 4 && choose4 < 2 &&  //
           choose1 >= 0 && choose2 >= 0 && choose3 >= 0 && choose4 >= 0);

    double dsignR  = (R[0] < 0) ? -1.0e0 : 1.0e0;
    double expB1R  = (R[0] < 0) ? exp(B1 * R[0]) : exp(-B1 * R[0]);
    double expB2RR = exp(-B2 * R[0] * R[0]);
    double expB3RR = exp(-B3L[choose4] * R[0] * R[0]);

    if (choose1 == 0) {
        V[0]   = E1;
        dV[0]  = 0.0e0;
        ddV[0] = 0.0e0;
    }
    if (choose1 == 1) {
        V[0]   = E2 + dsignR * A1 * (1 - expB1R);
        dV[0]  = A1 * B1 * expB1R;
        ddV[0] = -A1 * B1 * B1 * dsignR * expB1R;
    }
    if (choose1 == 2) {
        V[0]   = -E2 - dsignR * A1 * (1 - expB1R);
        dV[0]  = -A1 * B1 * expB1R;
        ddV[0] = A1 * dsignR * B1 * B1 * expB1R;
    }
    if (choose1 == 3) {
        V[0]   = E1 + A2 * expB2RR;
        dV[0]  = -2 * B2 * R[0] * expB2RR;
        ddV[0] = (2 * B2 * R[0] * R[0] - 1) * 2 * B2 * expB2RR;
    }
    if (choose1 == 4) {
        V[0]   = E1 - A2 * expB2RR;
        dV[0]  = 2 * B2 * R[0] * expB2RR;
        ddV[0] = -(2 * B2 * R[0] * R[0] - 1) * 2 * B2 * expB2RR;
    }

    if (choose2 == 0) {
        V[3]   = E2;
        dV[3]  = 0.0e0;
        ddV[3] = 0.0e0;
    }
    if (choose2 == 1) {
        V[3]   = E3 + dsignR * A1 * (1 - expB1R);
        dV[3]  = A1 * B1 * expB1R;
        ddV[3] = -A1 * B1 * B1 * dsignR * expB1R;
    }
    if (choose2 == 2) {
        V[3]   = E1 - dsignR * A1 * (1 - expB1R);
        dV[3]  = -A1 * B1 * expB1R;
        ddV[3] = A1 * dsignR * B1 * B1 * expB1R;
    }
    if (choose2 == 3) {
        V[3]   = E2 + A2 * expB2RR;
        dV[3]  = -2 * B2 * R[0] * expB2RR;
        ddV[3] = (2 * B2 * R[0] * R[0] - 1) * 2 * B2 * expB2RR;
    }
    if (choose2 == 4) {
        V[3]   = E2 - A2 * expB2RR;
        dV[3]  = 2 * B2 * R[0] * expB2RR;
        ddV[3] = -(2 * B2 * R[0] * R[0] - 1) * 2 * B2 * expB2RR;
    }

    V[1]   = A3L[choose3] * expB3RR;
    dV[1]  = -2 * B3L[choose4] * R[0] * expB3RR;
    ddV[1] = (2 * B3L[choose4] * R[0] * R[0] - 1) * 2 * B3L[choose4] * expB3RR;
    V[2]   = V[1];
    dV[2]  = dV[1];
    ddV[2] = ddV[1];
    return 0;
}

/**
 * @brief      MORSE3A model with parameter set-A
 * doi:10.1016/S0009-2614(01)01242-8
 */
int mspes_MORSE3A(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double De[3]    = {0.003e0, 0.004e0, 0.003e0};
    const double Be[3]    = {0.65e0, 0.60e0, 0.65e0};
    const double Re[3]    = {5.0e0, 4.0e0, 6.0e0};
    const double C[3]     = {0.0e0, 0.01e0, 0.006e0};
    const double Aijx[3]  = {0.002e0, 0.00e0, 0.002e0};
    const double Aeijx[3] = {16.0e0, 0.00e0, 16.0e0};
    const double Rijx[3]  = {4.80e0, 0.00e0, 3.40e0};

    double exp_bR1 = exp(-Be[0] * (R[0] - Re[0]));
    double exp_bR2 = exp(-Be[1] * (R[0] - Re[1]));
    double exp_bR3 = exp(-Be[2] * (R[0] - Re[2]));

    V[0] = De[0] * (1.0e0 - exp_bR1) * (1.0e0 - exp_bR1) + C[0];
    V[4] = De[1] * (1.0e0 - exp_bR2) * (1.0e0 - exp_bR2) + C[1];
    V[8] = De[2] * (1.0e0 - exp_bR3) * (1.0e0 - exp_bR3) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0e0 - exp_bR1) * exp_bR1 * 2 * Be[0];
    dV[4] = De[1] * (1.0e0 - exp_bR2) * exp_bR2 * 2 * Be[1];
    dV[8] = De[2] * (1.0e0 - exp_bR3) * exp_bR3 * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    ddV[0] = De[0] * (2 * exp_bR1 - 1.0e0) * exp_bR1 * 4 * Be[0] * Be[0];
    ddV[4] = De[1] * (2 * exp_bR2 - 1.0e0) * exp_bR2 * 4 * Be[1] * Be[1];
    ddV[8] = De[2] * (2 * exp_bR3 - 1.0e0) * exp_bR3 * 4 * Be[2] * Be[2];

    ddV[1] = dV[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2])) + V[1] * (-2 * Aeijx[2]);  // ddV(0,1)
    ddV[2] = dV[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1])) + V[2] * (-2 * Aeijx[1]);  // ddV(0,2)
    ddV[5] = dV[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0])) + V[5] * (-2 * Aeijx[0]);  // ddV(1,2)
    ddV[3] = ddV[1];                                                               // ddV(1,0)
    ddV[6] = ddV[2];                                                               // ddV(2,0)
    ddV[7] = ddV[5];                                                               // ddV(2,1)
    return 0;
}

/**
 * @brief      MORSE3B model with parameter set-B
 * doi:10.1016/S0009-2614(01)01242-8
 */
int mspes_MORSE3B(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double De[3]    = {0.020e0, 0.010e0, 0.003e0};
    const double Be[3]    = {0.65e0, 0.40e0, 0.65e0};
    const double Re[3]    = {4.5e0, 4.0e0, 4.4e0};
    const double C[3]     = {0.0e0, 0.01e0, 0.02e0};
    const double Aijx[3]  = {0.00e0, 0.005e0, 0.005e0};
    const double Aeijx[3] = {0.00e0, 32.0e0, 32.0e0};
    const double Rijx[3]  = {0.00e0, 3.34e0, 3.66e0};

    double exp_bR1 = exp(-Be[0] * (R[0] - Re[0]));
    double exp_bR2 = exp(-Be[1] * (R[0] - Re[1]));
    double exp_bR3 = exp(-Be[2] * (R[0] - Re[2]));

    V[0] = De[0] * (1.0e0 - exp_bR1) * (1.0e0 - exp_bR1) + C[0];
    V[4] = De[1] * (1.0e0 - exp_bR2) * (1.0e0 - exp_bR2) + C[1];
    V[8] = De[2] * (1.0e0 - exp_bR3) * (1.0e0 - exp_bR3) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0e0 - exp_bR1) * exp_bR1 * 2 * Be[0];
    dV[4] = De[1] * (1.0e0 - exp_bR2) * exp_bR2 * 2 * Be[1];
    dV[8] = De[2] * (1.0e0 - exp_bR3) * exp_bR3 * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    ddV[0] = De[0] * (2 * exp_bR1 - 1.0e0) * exp_bR1 * 4 * Be[0] * Be[0];
    ddV[4] = De[1] * (2 * exp_bR2 - 1.0e0) * exp_bR2 * 4 * Be[1] * Be[1];
    ddV[8] = De[2] * (2 * exp_bR3 - 1.0e0) * exp_bR3 * 4 * Be[2] * Be[2];

    ddV[1] = dV[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2])) + V[1] * (-2 * Aeijx[2]);  // ddV(0,1)
    ddV[2] = dV[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1])) + V[2] * (-2 * Aeijx[1]);  // ddV(0,2)
    ddV[5] = dV[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0])) + V[5] * (-2 * Aeijx[0]);  // ddV(1,2)
    ddV[3] = ddV[1];                                                               // ddV(1,0)
    ddV[6] = ddV[2];                                                               // ddV(2,0)
    ddV[7] = ddV[5];                                                               // ddV(2,1)
    return 0;
}

/**
 * @brief      MORSE3C model with parameter set-C
 * doi:10.1016/S0009-2614(01)01242-8
 */
int mspes_MORSE3C(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double De[3]    = {0.020e0, 0.020e0, 0.003e0};
    const double Be[3]    = {0.40e0, 0.65e0, 0.65e0};
    const double Re[3]    = {4.0e0, 4.5e0, 6.0e0};
    const double C[3]     = {0.02e0, 0.00e0, 0.02e0};
    const double Aijx[3]  = {0.00e0, 0.005e0, 0.005e0};
    const double Aeijx[3] = {0.00e0, 32.0e0, 32.0e0};
    const double Rijx[3]  = {0.00e0, 4.97e0, 3.40e0};

    double exp_bR1 = exp(-Be[0] * (R[0] - Re[0]));
    double exp_bR2 = exp(-Be[1] * (R[0] - Re[1]));
    double exp_bR3 = exp(-Be[2] * (R[0] - Re[2]));

    V[0] = De[0] * (1.0e0 - exp_bR1) * (1.0e0 - exp_bR1) + C[0];
    V[4] = De[1] * (1.0e0 - exp_bR2) * (1.0e0 - exp_bR2) + C[1];
    V[8] = De[2] * (1.0e0 - exp_bR3) * (1.0e0 - exp_bR3) + C[2];

    V[1] = Aijx[2] * exp(-Aeijx[2] * (R[0] - Rijx[2]) * (R[0] - Rijx[2]));  // V(0,1)
    V[2] = Aijx[1] * exp(-Aeijx[1] * (R[0] - Rijx[1]) * (R[0] - Rijx[1]));  // V(0,2)
    V[5] = Aijx[0] * exp(-Aeijx[0] * (R[0] - Rijx[0]) * (R[0] - Rijx[0]));  // V(1,2)
    V[3] = V[1];                                                            // V(1,0)
    V[6] = V[2];                                                            // V(2,0)
    V[7] = V[5];                                                            // V(2,1)

    if (flag < 1) return 0;

    dV[0] = De[0] * (1.0e0 - exp_bR1) * exp_bR1 * 2 * Be[0];
    dV[4] = De[1] * (1.0e0 - exp_bR2) * exp_bR2 * 2 * Be[1];
    dV[8] = De[2] * (1.0e0 - exp_bR3) * exp_bR3 * 2 * Be[2];

    dV[1] = V[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2]));  // dV(0,1)
    dV[2] = V[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1]));  // dV(0,2)
    dV[5] = V[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0]));  // dV(1,2)
    dV[3] = dV[1];                                      // dV(1,0)
    dV[6] = dV[2];                                      // dV(2,0)
    dV[7] = dV[5];                                      // dV(2,1)

    if (flag < 2) return 0;
    ddV[0] = De[0] * (2 * exp_bR1 - 1.0e0) * exp_bR1 * 4 * Be[0] * Be[0];
    ddV[4] = De[1] * (2 * exp_bR2 - 1.0e0) * exp_bR2 * 4 * Be[1] * Be[1];
    ddV[8] = De[2] * (2 * exp_bR3 - 1.0e0) * exp_bR3 * 4 * Be[2] * Be[2];

    ddV[1] = dV[1] * (-2 * Aeijx[2] * (R[0] - Rijx[2])) + V[1] * (-2 * Aeijx[2]);  // ddV(0,1)
    ddV[2] = dV[2] * (-2 * Aeijx[1] * (R[0] - Rijx[1])) + V[2] * (-2 * Aeijx[1]);  // ddV(0,2)
    ddV[5] = dV[5] * (-2 * Aeijx[0] * (R[0] - Rijx[0])) + V[5] * (-2 * Aeijx[0]);  // ddV(1,2)
    ddV[3] = ddV[1];                                                               // ddV(1,0)
    ddV[6] = ddV[2];                                                               // ddV(2,0)
    ddV[7] = ddV[5];                                                               // ddV(2,1)
    return 0;
}

/**
 * @brief      MORSE15 model
 * doi:10.1063/1.4742155
 */
int mspes_MORSE15(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double Dg = 0.2e0, De = 0.05e0, Dc = 0.05e0, alpha = 0.4e0, lambda = 1.0e0, eta = 0.004e0;
    const double Re[15] = {3.0000000000000e0,    6.8335793401926175e0, 6.909940519903887e0, 6.988372712202933e0,
                           7.0690037103879675e0, 7.151973244334301e0,  7.237434522736795e0, 7.325556032471846e0,
                           7.416523648311893e0,  7.5105431195440095e0, 7.607843017311827e0, 7.708678249086644e0,
                           7.813334276495011e0,  7.922132212503321e0,  8.035435027584366};
    for (int i = 0; i < 225; ++i) V[i] = 0.0e0;
    V[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * (1 - exp(-alpha * (R[0] - Re[0])));
    for (int i = 1; i < 15; ++i) {
        V[i * 16] = exp(-alpha * R[0]) + eta * (i + 1) + De;
        V[i]      = Dc * exp(-lambda * (R[0] - Re[i]) * (R[0] - Re[i]));
        V[i * 15] = V[i];
    }
    if (flag < 1) return 0;

    for (int i = 0; i < 225; ++i) dV[i] = 0.0e0;
    dV[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * exp(-alpha * (R[0] - Re[0])) * 2 * alpha;
    for (int i = 1; i < 15; ++i) {
        dV[i * 16] = -alpha * exp(-alpha * R[0]);
        dV[i]      = V[i] * (-2 * lambda * (R[0] - Re[i]));
        dV[i * 15] = dV[i];
    }

    if (flag < 2) return 0;
    for (int i = 0; i < 225; ++i) ddV[i] = 0.0e0;
    ddV[0] = Dg * (2.0 * exp(-alpha * (R[0] - Re[0])) - 1) * exp(-alpha * (R[0] - Re[0])) * 4 * alpha * alpha;
    for (int i = 1; i < 15; ++i) {
        ddV[i * 16] = alpha * alpha * exp(-alpha * R[0]);
        ddV[i]      = dV[i] * (-2 * lambda * (R[0] - Re[i])) + V[i] * (-2 * lambda);
        ddV[i * 15] = ddV[i];
    }
    return 0;
}

/**
 * @brief      MORSE15 model revised with constant coupling
 */
int mspes_MORSE15C(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double Dg = 0.2e0, De = 0.05e0, Dc = 0.05e0, alpha = 0.4e0, lambda = 1.0e0, eta = 0.004e0;
    const double Re[15] = {3.0000000000000e0,    6.8335793401926175e0, 6.909940519903887e0, 6.988372712202933e0,
                           7.0690037103879675e0, 7.151973244334301e0,  7.237434522736795e0, 7.325556032471846e0,
                           7.416523648311893e0,  7.5105431195440095e0, 7.607843017311827e0, 7.708678249086644e0,
                           7.813334276495011e0,  7.922132212503321e0,  8.035435027584366};
    const double kappa  = 0.25e0;

    for (int i = 0; i < 225; ++i) V[i] = 0.0e0;
    V[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * (1 - exp(-alpha * (R[0] - Re[0])));
    for (int i = 1; i < 15; ++i) {
        V[i * 16] = exp(-alpha * R[0]) + eta * (i + 1) + De;
        V[i]      = kappa * Dc;
        V[i * 15] = V[i];
    }
    if (flag < 1) return 0;

    for (int i = 0; i < 225; ++i) dV[i] = 0.0e0;
    dV[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * exp(-alpha * (R[0] - Re[0])) * 2 * alpha;
    for (int i = 1; i < 15; ++i) {
        dV[i * 16] = -alpha * exp(-alpha * R[0]);
        // dV[i]      = 0.0e0;
        // dV[i * 15] = dV[i];
    }

    if (flag < 2) return 0;
    for (int i = 0; i < 225; ++i) ddV[i] = 0.0e0;
    ddV[0] = Dg * (2.0 * exp(-alpha * (R[0] - Re[0])) - 1) * exp(-alpha * (R[0] - Re[0])) * 4 * alpha * alpha;
    for (int i = 1; i < 15; ++i) {
        ddV[i * 16] = alpha * alpha * exp(-alpha * R[0]);
        // ddV[i]      = 0.0e0;
        // ddV[i * 15] = ddV[i];
    }
    return 0;
}

/**
 * @brief      MORSE15 model revised with exponential coupling
 */
int mspes_MORSE15E(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    // V
    const double Dg = 0.2e0, De = 0.05e0, Dc = 0.05e0, alpha = 0.4e0, lambda = 1.0e0, eta = 0.004e0;
    const double Re[15] = {3.0000000000000e0,    6.8335793401926175e0, 6.909940519903887e0, 6.988372712202933e0,
                           7.0690037103879675e0, 7.151973244334301e0,  7.237434522736795e0, 7.325556032471846e0,
                           7.416523648311893e0,  7.5105431195440095e0, 7.607843017311827e0, 7.708678249086644e0,
                           7.813334276495011e0,  7.922132212503321e0,  8.035435027584366};
    const double kappa  = 0.2e0;

    for (int i = 0; i < 225; ++i) V[i] = 0.0e0;
    V[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * (1 - exp(-alpha * (R[0] - Re[0])));
    for (int i = 1; i < 15; ++i) {
        V[i * 16] = exp(-alpha * R[0]) + eta * (i + 1) + De;
        V[i]      = Dc * exp(-kappa * R[0]);
        V[i * 15] = V[i];
    }
    if (flag < 1) return 0;

    for (int i = 0; i < 225; ++i) dV[i] = 0.0e0;
    dV[0] = Dg * (1 - exp(-alpha * (R[0] - Re[0]))) * exp(-alpha * (R[0] - Re[0])) * 2 * alpha;
    for (int i = 1; i < 15; ++i) {
        dV[i * 16] = -alpha * exp(-alpha * R[0]);
        dV[i]      = -kappa * V[i];
        dV[i * 15] = dV[i];
    }

    if (flag < 2) return 0;
    for (int i = 0; i < 225; ++i) ddV[i] = 0.0e0;
    ddV[0] = Dg * (2.0 * exp(-alpha * (R[0] - Re[0])) - 1) * exp(-alpha * (R[0] - Re[0])) * 4 * alpha * alpha;
    for (int i = 1; i < 15; ++i) {
        ddV[i * 16] = alpha * alpha * exp(-alpha * R[0]);
        ddV[i]      = kappa * kappa * V[i];
        ddV[i * 15] = ddV[i];
    }
    return 0;
}

/**
 * @brief      Caldeira-Leggett model (single-mode spin-boson model)
 * doi:
 * html:http://www.scholarpedia.org/article/Caldeira-Leggett_model
 */
int mspes_CL1D(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    double& bias  = mspes_parm_real[0];
    double& delta = mspes_parm_real[1];
    double& coup  = mspes_parm_real[2];
    double& wmod  = mspes_parm_real[3];

    double meanV = 0.5e0 * wmod * wmod * R[0] * R[0];  // mass = 1
    V[0]         = meanV + bias + coup * R[0];
    V[3]         = meanV - bias - coup * R[0];
    V[1]         = delta;
    V[2]         = V[1];
    if (flag < 1) return 0;

    dV[0] = wmod * wmod * R[0] + coup;
    dV[3] = wmod * wmod * R[0] - coup;
    dV[1] = 0.0e0;
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = wmod * wmod;
    ddV[3] = ddV[0];
    dV[1]  = 0.0e0;
    dV[2]  = dV[1];
    return 0;
}

/**
 * @brief      Jaynes–Cummings model
 * doi:
 */
int mspes_JC1D(double* V, double* dV, double* ddV, double* R, double* P, int flag, int rdim, int fdim) {
    double& bias  = mspes_parm_real[0];
    double& delta = mspes_parm_real[1];
    double& coup  = mspes_parm_real[2];
    double& wmod  = mspes_parm_real[3];
    double  sqrtw = sqrt(wmod);

    double meanV = 0.5e0 * wmod * wmod * R[0] * R[0];  // mass = 1
    V[0]         = meanV + bias;
    V[3]         = meanV - bias;
    V[1]         = delta + coup * (sqrtw * R[0]
                           //+ phys::math::im * P[0] / sqrtw
                          );
    V[2]         = V[1];  // std::conj(V[1]);
    if (flag < 1) return 0;

    dV[0] = wmod * wmod * R[0];
    dV[3] = wmod * wmod * R[0];
    dV[1] = coup;
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = wmod * wmod;
    ddV[3] = ddV[0];
    dV[1]  = 0.0e0;
    dV[2]  = dV[1];
    return 0;
}

/**
 * @brief      Rabi model (transformed gauge)
 * doi:
 */
int mspes_RABI(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    double& bias   = mspes_parm_real[0];
    double& delta  = mspes_parm_real[1];
    double& coup   = mspes_parm_real[2];
    double& wmod   = mspes_parm_real[3];
    double  sqrt2w = sqrt(2.0e0 * wmod);

    double meanV = 0.5e0 * wmod * wmod * R[0] * R[0];  // mass = 1
    V[0]         = meanV + bias;
    V[3]         = meanV - bias;
    V[1]         = delta + sqrt2w * coup * R[0];
    V[2]         = V[1];
    if (flag < 1) return 0;

    dV[0] = wmod * wmod * R[0];
    dV[3] = wmod * wmod * R[0];
    dV[1] = sqrt2w * coup;
    dV[2] = dV[1];
    if (flag < 2) return 0;

    // ddV
    ddV[0] = wmod * wmod;
    ddV[3] = ddV[0];
    dV[1]  = 0.0e0;
    dV[2]  = dV[1];
    return 0;
}

int mspes_NA_I(double* V, double* dV, double* ddV, double* R, int flag, int rdim, int fdim) {
    const double au_2_ev  = 27.21138602e0;
    const double au_2_ang = 0.529177249e0;
    const double Acov     = 3150.0e0 / au_2_ev;
    const double Bcov12   = pow(2.647e0 / au_2_ang, 12) / au_2_ev;
    const double Ccov     = 1000.0e0 / au_2_ev / pow(au_2_ang, 6);
    const double Rcov     = 0.435e0 / au_2_ang;
    const double Aion     = 2760.0e0 / au_2_ev;
    const double Bion8    = pow(2.398e0 / au_2_ang, 8) / au_2_ev;
    const double Cion     = 11.3e0 / au_2_ev / pow(au_2_ang, 6);
    const double Rion     = 0.3489e0 / au_2_ang;
    const double al_Mp    = 0.408e0 / pow(au_2_ang, 3);
    const double al_Xm    = 6.431e0 / pow(au_2_ang, 3);
    const double al_M     = 27.0e0 / pow(au_2_ang, 3);
    const double al_X     = 7.0e0 / pow(au_2_ang, 3);
    const double Eth      = 2.075e0 / au_2_ev;
    const double Aoff     = 17.08e0 / au_2_ev;
    const double Roff     = 1.239e0 / au_2_ang;
    const double mred     = (22.989769282e0 * 126.904473e0) / (22.989769282e0 + 126.904473e0) / phys::au_2_amu;

    double r   = R[0];
    double r2  = r * r;
    double r3  = r2 * r;
    double r4  = r2 * r2;
    double r5  = r4 * r;
    double r6  = r2 * r4;
    double r7  = r6 * r;
    double r8  = r4 * r4;
    double r9  = r8 * r;
    double r10 = r9 * r;
    double r12 = r6 * r6;
    double r13 = r12 * r;
    double r14 = r13 * r;

    double Vcent = mspes_parm_real[1] / (2 * mred * r2);

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

    double dVcent = -2 * mspes_parm_real[1] / (2 * mred * r2 * r);

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

    double ddVcent = 6 * mspes_parm_real[1] / (2 * mred * r4);

    ddV[0] = ddVcent                                                       //
             + (156 * Bcov12 / r14) * exp(-r / Rcov)                       //
             + 2 * (-12 * Bcov12 / r13) * (-1 / Rcov) * exp(-r / Rcov)     //
             + (Acov + Bcov12 / r12) * 1 / (Rcov * Rcov) * exp(-r / Rcov)  //
             - 42 * Ccov / r8;                                             //
    ddV[3] = ddVcent                                                       //
             + (72 * Bion8 / r10) * exp(-r / Rion)                         //
             + 2 * (-8 * Bion8 / r9) * (-1 / Rion) * exp(-r / Rion)        //
             + (Aion + Bion8 / r8) * 1 / (Rion * Rion) * exp(-r / Rion)    //
             - 2 / r3                                                      //
             - 10 * (al_Mp + al_Xm) / r6                                   //
             - 42 * Cion / r8                                              //
             - 112 * al_Mp * al_Xm / r9;                                   //
    ddV[1] = 1 / (Roff * Roff) * Aoff * exp(-r / Roff);
    ddV[2] = ddV[1];
    return 0;
}

const std::string Model_NAD1D::getName() { return "Model_NAD1D"; }

int Model_NAD1D::getType() const { return utils::hash(FUNCTION_NAME); }

void Model_NAD1D::setInputParam_impl(std::shared_ptr<Param> PM) {
    nad1d_type = NAD1DPolicy::_from(_param->get_string({"model.nad1d_flag"}, LOC(), "SAC"));

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
            double x0_read   = _param->get_real({"model.x0"}, LOC(), -10.0e0);
            double p0_read   = _param->get_real({"model.p0"}, LOC(), 10.0e0);
            double mass_read = _param->get_real({"model.m0"}, LOC(), 2000.0e0);
            double tend_read = _param->get_real({"model.tend"}, LOC(), -1);
            double dt_read   = _param->get_real({"model.dt"}, LOC(), -1);

            double tend_rev = std::abs(5 * x0_read * mass_read / p0_read);
            double dt_rev   = std::abs(x0_read * mass_read / p0_read) / 5000;
            if (tend_read < 0) _param->set_real("model.tend", tend_rev);
            if (dt_read < 0) _param->set_real("model.dt", dt_rev);
            break;
        }
        case NAD1DPolicy::NA_I: {
            // CHECK_EQ(F, 2);
            // initializetion in au
            double mass_Na = 22.989769282e0 / phys::au_2_amu;
            double mass_I  = 126.904473e0 / phys::au_2_amu;
            double mred    = (mass_Na * mass_I) / (mass_Na + mass_I);

            double parm_E    = _param->get_real({"model.parm_E"}, LOC(), 1.0e0);
            double tend_read = _param->get_real({"model.tend"}, LOC(), -1);
            double dt_read   = _param->get_real({"model.dt"}, LOC(), -1);

            double x0_max   = 160.0e0 / phys::au_2_ang;
            double vel      = sqrt(2 * parm_E / mred);
            double tend_rev = 2 * x0_max / vel;
            double dt_rev   = tend_rev / 50000;
            if (tend_read < 0) _param->set_real("model.tend", tend_rev);
            if (dt_read < 0) _param->set_real("model.dt", dt_rev);
            break;
        }
    }
};

void Model_NAD1D::setInputDataSet_impl(std::shared_ptr<DataSet> DS) {
    Hsys = DS->def(DATA::model::Hsys);
    memset(Hsys.data(), 0, Dimension::FF * sizeof(kids_real));

    if (nad1d_type == NAD1DPolicy::PURE) {
        std::ifstream ifs("Hsys.dat");
        std::string   H_unit_str;
        std::string   firstline;
        getline(ifs, firstline);
        std::stringstream sstr(firstline);
        sstr >> H_unit_str;
        double    H_unit = phys::us::conv(phys::au::unit, phys::us::parse(H_unit_str));
        kids_real val;
        for (int i = 0; i < Dimension::FF; ++i)
            if (ifs >> val) Hsys[i] = val / H_unit;
        ifs.close();
    }

    // model field
    mass = DS->def(DATA::model::mass);
    vpes = DS->def(DATA::model::vpes);
    grad = DS->def(DATA::model::grad);
    hess = DS->def(DATA::model::hess);
    V    = DS->def(DATA::model::V);
    dV   = DS->def(DATA::model::dV);
    // ddV  = DS->def(DATA::model::ddV);

    x0      = DS->def(DATA::model::x0);
    p0      = DS->def(DATA::model::p0);
    x_sigma = DS->def(DATA::model::x_sigma);
    p_sigma = DS->def(DATA::model::p_sigma);

    // init & integrator
    x      = DS->def(DATA::integrator::x);
    p      = DS->def(DATA::integrator::p);
    p_sign = DS->def(DATA::integrator::p_sign);

    double  x0_read = _param->get_real({"model.x0grid"}, LOC(), -10.0e0);
    int     Nxgird  = _param->get_int({"model.Nxgrid"}, LOC(), 101);
    double  dx      = (2 * fabs(x0_read)) / (Nxgird - 1);
    double* xgrid   = DS->def_real("integrator.xgrid", Nxgird);
    for (int i = 0; i < Nxgird; ++i) xgrid[i] = -abs(x0_read) + i * dx;

    mass[0]     = _param->get_real({"model.m0"}, LOC(), 2000.0e0);
    x0[0]       = _param->get_real({"model.x0"}, LOC(), 100.0e0);
    p0[0]       = _param->get_real({"model.p0"}, LOC(), 100.0e0);
    double varx = _param->get_real({"model.varx"}, LOC(), 0.5e0);
    double varp = _param->get_real({"model.varp"}, LOC(), 0.5e0);
    x_sigma[0]  = sqrt(varx);
    p_sigma[0]  = sqrt(varp);

    switch (nad1d_type) {
        case NAD1DPolicy::SAC:
        case NAD1DPolicy::SAC2:
        case NAD1DPolicy::DAC:
        case NAD1DPolicy::ECR:
        case NAD1DPolicy::DBG:
        case NAD1DPolicy::DAG:
        case NAD1DPolicy::DRN: {
            mass[0] = 2000.0e0;
            break;
        }
        case NAD1DPolicy::DPES: {
            mass[0]           = 2000.0e0;
            mspes_parm_int[0] = _param->get_real({"model.dpes.choose1"}, LOC(), 1);
            mspes_parm_int[1] = _param->get_real({"model.dpes.choose2"}, LOC(), 0);
            mspes_parm_int[2] = _param->get_real({"model.dpes.choose3"}, LOC(), 1);
            mspes_parm_int[3] = _param->get_real({"model.dpes.choose4"}, LOC(), 0);
            break;
        }
        case NAD1DPolicy::SAC3: {  // asymmetrical SAC
            mass[0]           = 1980.0e0;
            double gammawidth = 0.25e0;
            x_sigma[0]        = sqrt(0.5e0 / gammawidth);
            p_sigma[0]        = sqrt(0.5e0 * gammawidth);
            break;
        }
        case NAD1DPolicy::MORSE3A:
        case NAD1DPolicy::MORSE3B:
        case NAD1DPolicy::MORSE3C: {
            // CHECK_EQ(F, 3);
            mass[0] = 20000.0e0;
            switch (nad1d_type) {
                case NAD1DPolicy::MORSE3A:
                    x0[0] = 2.9e0;
                    break;
                case NAD1DPolicy::MORSE3B:
                    x0[0] = 3.3e0;
                    break;
                case NAD1DPolicy::MORSE3C:
                    x0[0] = 2.1e0;
                    break;
            }
            p0[0]          = 0.0e0;
            double wground = 5.0e-03;
            x_sigma[0]     = sqrt(0.5e0 / (mass[0] * wground));
            p_sigma[0]     = 0.5e0 / x_sigma[0];
            break;
        }
        case NAD1DPolicy::MORSE15:
        case NAD1DPolicy::MORSE15C:
        case NAD1DPolicy::MORSE15E: {
            // CHECK_EQ(F, 15);
            mass[0]   = 2000.0e0;
            x0[0]     = 13.0e0;
            p0[0]     = -30.0e0;
            double Dg = 0.2e0, alpha = 0.4e0;  // 2*alpha^2 * D = m*w^2
            x_sigma[0] = sqrt(0.5e0 / std::sqrt(mass[0] * 2 * alpha * alpha * Dg));
            p_sigma[0] = 0.5e0 / x_sigma[0];
            break;
        }
        case NAD1DPolicy::CL1D:
        case NAD1DPolicy::JC1D: {
            // CHECK_EQ(F, 2);
            mass[0]            = 1.0e0;
            mspes_parm_real[0] = _param->get_real({"model.nad1d_e"}, LOC(), 1.0e0);
            mspes_parm_real[1] = _param->get_real({"model.nad1d_d"}, LOC(), 1.0e0);
            mspes_parm_real[2] = _param->get_real({"model.nad1d_c"}, LOC(), 1.0e0);
            mspes_parm_real[3] = _param->get_real({"model.nad1d_w"}, LOC(), 1.0e0);
            break;
        }
        case NAD1DPolicy::NA_I: {
            // CHECK_EQ(F, 2);
            // initializetion in au
            double mass_Na = 22.989769282e0 / phys::au_2_amu;
            double mass_I  = 126.904473e0 / phys::au_2_amu;
            mass[0]        = (mass_Na * mass_I) / (mass_Na + mass_I);

            mspes_parm_real[0] = _param->get_real({"model.parm_E"}, LOC(), 1.0e0);
            x0[0]              = 160.0e0 / phys::au_2_ang;
            p0[0]              = -sqrt(2 * mass[0] * mspes_parm_real[0]);
            break;
        }
    }
};

Status& Model_NAD1D::initializeKernel_impl(Status& stat) {
    // executeKernel(stat);
    return stat;  // @todo

    switch (nad1d_type) {
        case NAD1DPolicy::NA_I: {
            double bmax = 9.0e0 / phys::au_2_ang;
            double randu;
            Kernel_Random::rand_uniform(&randu);
            double b           = sqrt(randu) * bmax;
            mspes_parm_real[1] = 2 * mass[0] * mspes_parm_real[0] * b * b;  // l^2 = 2 m E b^2
            for (int j = 0; j < Dimension::N; ++j) {
                x[j] = x0[j];
                p[j] = p0[j];
            }
            break;
        }
        default: {
            Kernel_Random::rand_gaussian(x.data(), Dimension::N);
            Kernel_Random::rand_gaussian(p.data(), Dimension::N);
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
    _dataset->def(DATA::init::x, x);
    _dataset->def(DATA::init::p, p);
    executeKernel(stat);
    return stat;
}

Status& Model_NAD1D::executeKernel_impl(Status& stat) {
    switch (nad1d_type) {
        case NAD1DPolicy::SAC:
            mspes_SAC(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::SAC2:
            mspes_SAC2(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::SAC3:
            mspes_SAC3(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::DAC:
            mspes_DAC(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::ECR:
            mspes_ECR(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::DBG:
            mspes_DBG(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::DAG:
            mspes_DAG(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::DRN:
            mspes_DRN(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::DPES:
            mspes_DPES(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3A:
            mspes_MORSE3A(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3B:
            mspes_MORSE3B(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE3C:
            mspes_MORSE3C(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE15:
            mspes_MORSE15(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE15C:
            mspes_MORSE15C(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::MORSE15E:
            mspes_MORSE15E(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::JC1D:
            mspes_JC1D(V.data(), dV.data(), ddV.data(), x.data(), p.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::RABI:
            mspes_RABI(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::CL1D:
            mspes_CL1D(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::NA_I:
            mspes_NA_I(V.data(), dV.data(), ddV.data(), x.data(), 1, Dimension::N, Dimension::F);
            break;
        case NAD1DPolicy::PURE: {
            if (count_exec == 0) {
                for (int i = 0; i < Dimension::FF; ++i) { V[i] = Hsys[i], dV[i] = 0.0e0; }
            }
            break;
        }
    }
    if (p[0] >= 0) {
        p_sign[0] = phys::math::iu, p_sign[1] = phys::math::iz;
    } else {
        p_sign[0] = phys::math::iz, p_sign[1] = phys::math::iu;
    }
    return stat;
}


};  // namespace PROJECT_NS
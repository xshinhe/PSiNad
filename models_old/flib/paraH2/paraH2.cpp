#include <cmath>

// Liquid H2 (Ref https://aip.scitation.org/doi/pdf/10.1063/1.437103)

int paraH2_pair(double *V, double *dV, double *r) {
    constexpr double a0  = 1.713f;
    constexpr double a1  = 1.5671f;
    constexpr double a2  = 0.00993f;
    constexpr double c6  = 12.14f;
    constexpr double c8  = 215.2f;
    constexpr double c9  = 143.1f;
    constexpr double c10 = 4813.9f;

    constexpr double rc = 8.32;  // = 1.28*rm (rm=3.41A), ?@tocheck


    double r2  = r[0] * r[0];
    double r6  = r2 * r2 * r2;
    double r8  = r6 * r2;
    double r9  = r8 * r[0];
    double r10 = r8 * r2;

    double vpes_exp = exp(a0 - a1 * r[0] - a2 * r2);
    double vpes_pow = (-c6 / r6 - c8 / r8 + c9 / r9 - c10 / r10);

    if (r[0] < rc) {
        double z  = rc / r[0] - 1;
        double fc = exp(-z * z);

        V[0]  = vpes_exp + vpes_pow * fc;
        dV[0] = (-a1 - 2 * a2 * r[0]) * vpes_exp                                          // term1
                + (6 * c6 / r6 + 8 * c8 / r8 - 9 * c9 / r9 + 10 * c10 / r10) / r[0] * fc  // term2-a
                + vpes_pow * fc * (-2 * z) * (-rc / r2);                                  // term2-b
    } else {
        V[0]  = vpes_exp + vpes_pow;
        dV[0] = (-a1 - 2 * a2 * r[0]) * vpes_exp                                      // term1
                + (6 * c6 / r6 + 8 * c8 / r8 - 9 * c9 / r9 + 10 * c10 / r10) / r[0];  // term2
    }
    return 0;
}

int paraH2_pot(double *V, double *dV, double *nr, const int &N, const double &BoxL) {
    V[0] = 0.0f;
    for (int i = 0; i < N; ++i) dV[i] = 0.0;

    // temperoray variables
    int Natom = N / 3;
    double r, dx[3];

    for (int iatom = 0, idxr1 = 0; iatom < Natom; ++iatom, idxr1 += 3) {
        double *rA1 = nr + idxr1;  // first atom
        double *fA1 = dV + idxr1;  // first atoms
        for (int jatom = iatom + 1, idxr2 = idxr1 + 3; jatom < Natom; ++jatom, idxr2 += 3) {
            double *rA2 = nr + idxr2;  // second atom
            double *fA2 = dV + idxr2;  // second atom

            r = 0.0f;
            for (int i = 0; i < 3; ++i) {
                dx[i] = rA1[i] - rA2[i];
                dx[i] -= round(dx[i] / BoxL) * BoxL;
                r += dx[i] * dx[i];
            }
            r = sqrt(r);

            constexpr double rcut = 16;
            if (r < rcut) {
                double Vij, dVij;
                paraH2_pair(&Vij, &dVij, &r);

                V[0] += Vij;
                for (int i = 0; i < 3; ++i) {
                    fA1[i] += dVij * dx[i] / r;
                    fA2[i] -= dVij * dx[i] / r;
                }
            }
        }
    }
    return 0;
}

#include "manysite_models.h"

using namespace ARRAY_EG;

ManySite_ForceField::ManySite_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    CHECK_EQ(N, 1);
    CHECK_EQ(F, 2);

    M     = Param_GetT(int, parm, "M", 100);
    Jp    = Param_GetT(double, parm, "Jp", 0.0f);
    Jz    = Param_GetT(double, parm, "Jz", 1.0f);
    alpha = Param_GetT(double, parm, "alpha", 0.0f);
    omega = Param_GetT(double, parm, "omega", 0.0f);

    /** Ising model and XY model
        Jij ≡ J*[1 − 3cos^2(θ)] / rij^α
    */
    ALLOCATE_PTR_TO_VECTOR(JpMat, M * M);
    ALLOCATE_PTR_TO_VECTOR(JzMat, M * M);
    ALLOCATE_PTR_TO_VECTOR(redX, M);
    ALLOCATE_PTR_TO_VECTOR(redY, M);
    ALLOCATE_PTR_TO_VECTOR(redZ, M);
    for (int is = 0, k = 0; is < M; ++is) {
        for (int js = 0; js < M; ++js, ++k) {
            if (is == js) {
                JpMat[k] = 0.0f, JzMat[k] = 0.0f;
                continue;
            }
            int dist = is > js ? is - js : js - is;
            JpMat[k] = Jp / pow(dist, alpha), JzMat[k] = Jz / pow(dist, alpha);
        }
    }
};

int ManySite_ForceField::ForceField_heff(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim) {
    ForceField_heff_Ising(H, rhos, mdim, fdim);
    return 0;
}

int ManySite_ForceField::ForceField_heff_Ising(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim) {
    memset(H, 0, M * FF * sizeof(num_complex));

    // reduce X/Y/Z in temporary array redX, redY, redZ (saving time)
    for (int is = 0, idxrho = 0; is < M; ++is, idxrho += FF) {
        num_complex* rhoi = rhos + idxrho;
        num_complex res;
        // X
        res = phys::math::iz;
        for (int i = 0, idx1 = 0; i < F; ++i) {
            for (int j = 0, idx2 = i; j < F; ++j, ++idx1, idx2 += F) { res += rhoi[idx1] * Pauli::X[idx2]; }
        }
        redX[is] = res;
        // Y
        res = phys::math::iz;
        for (int i = 0, idx1 = 0; i < F; ++i) {
            for (int j = 0, idx2 = i; j < F; ++j, ++idx1, idx2 += F) { res += rhoi[idx1] * Pauli::Y[idx2]; }
        }
        redY[is] = res;
        // Z
        res = phys::math::iz;
        for (int i = 0, idx1 = 0; i < F; ++i) {
            for (int j = 0, idx2 = i; j < F; ++j, ++idx1, idx2 += F) { res += rhoi[idx1] * Pauli::Z[idx2]; }
        }
        redZ[is] = res;
    }


    for (int is = 0, ij = 0; is < M; ++is) {
        num_complex* Hi = H + is * FF;

        // one-body interaction
        for (int k = 0; k < FF; ++k) Hi[k] = omega * Pauli::X[k];

        // two-body interaction
        for (int js = 0; js < M; ++js, ++ij) {
            for (int k = 0; k < FF; ++k) Hi[k] += 0.25e0 * Pauli::X[k] * JpMat[ij] * redX[js];
            for (int k = 0; k < FF; ++k) Hi[k] += 0.25e0 * Pauli::Y[k] * JpMat[ij] * redY[js];
            for (int k = 0; k < FF; ++k) Hi[k] += 0.5e0 * Pauli::Z[k] * JzMat[ij] * redZ[js];

            // three-body interaction: pass
            // ...
        }
    }
    return 0;
};

#include "forcefieldbase.h"

#include "../utils/definitions.h"

ForceField::ForceField(const Param& iparm) : Model(iparm){};
ForceField::~ForceField(){};

BO_ForceField::BO_ForceField(const Param& iparm) : ForceField(iparm) {
    // in parm (updated!) not in iparm
    Ndim     = Param_GetT(int, parm, "Ndim", 3);  // no. of dimension
    N        = Param_GetT(int, parm, "rdim");     // no. of nuclear DOFs
    NN       = N * N;
    mod_dt   = Param_GetT(double, parm, "dt", -1.0);
    mod_tend = Param_GetT(double, parm, "tend", -1.0);

    // CHECK_EQ(N % Ndim, 0);
    ALLOCATE_PTR_TO_VECTOR(mod_M, N);
    ALLOCATE_PTR_TO_VECTOR(mod_W, N);
    ALLOCATE_PTR_TO_VECTOR(mod_R0, N);
    ALLOCATE_PTR_TO_VECTOR(mod_P0, N);
    ALLOCATE_PTR_TO_VECTOR(mod_sigmaR, N);
    ALLOCATE_PTR_TO_VECTOR(mod_sigmaP, N);

    for (int i = 0; i < N; ++i) mod_R0[i] = 0.0f, mod_P0[i] = 0.0f, mod_M[i] = 1.0f, mod_W[i] = 0.0f;
}
BO_ForceField::~BO_ForceField(){};

int BO_ForceField::get_N() { return N; }
int BO_ForceField::get_Ndim() { return Ndim; }
double BO_ForceField::Suggest_dt() { return mod_dt; }
double BO_ForceField::Suggest_tend() { return mod_tend; }

int BO_ForceField::ForceField_init(double* nr, double* np, double* nm, const int& rdim, const int& itraj) { return 0; }

int BO_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim) {
    return (rdim > 0) ? 0 : 1;  // default, only 1 specification
}

int BO_ForceField::nspec() { return ForceField_spec(nullptr, nullptr, nullptr, -1); }

int BO_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                   const int& rdim) {
    V[0] = 0.0f;
    for (int i = 0; i < rdim; ++i) dV[i] = 0.0f;
    return 0;
}
int BO_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                   const int& rdim, const int& itraj, const int& isamp) {
    return ForceField_npes(V, dV, ddV, R, P, flag, rdim);  // default
}

int BO_ForceField::ForceField_init_default_build(const double& beta, const int& rdim) {
    double Qoverbeta;
    for (int j = 0; j < rdim; ++j) {
        /* note:
            for finite temperature: Qoverbeta = 0.5*W / dtanh(0.5*beta*W)
            for zero temperature: Qoverbeta = 0.5*W  // let beta very large or just be negative
        */
        Qoverbeta     = 0.5f * mod_W[j] / (beta > 0 ? std::tanh(0.5f * beta * mod_W[j]) : 1.0f);
        mod_sigmaR[j] = std::sqrt(Qoverbeta / (mod_M[j] * mod_W[j] * mod_W[j]));
        mod_sigmaP[j] = std::sqrt(Qoverbeta * mod_M[j]);
    }
    return 0;
};

// Harmonic Wigner/Wavepacket sampling.
int BO_ForceField::ForceField_init_default(double* nr, double* np, double* nm, const int& rdim, const int& itraj) {
    // LOG(INFO) << "calling bo default init"; @blame: log too many times
    rand_gaussian(nr, N), rand_gaussian(np, N);
    for (int i = 0; i < N; ++i) {
        nm[i] = mod_M[i];
        nr[i] = mod_R0[i] + nr[i] * mod_sigmaR[i];
        np[i] = mod_P0[i] + np[i] * mod_sigmaP[i];
    }
    return 0;
}

int BO_ForceField::CheckForceField() {
    LOG(INFO) << "After building " << std::endl << parm;
    double maxw = 0.0f, sumw = 0.0f;
    for (int i = 0; i < N; ++i) {
        if (std::abs(mod_W[i]) > maxw) maxw = mod_W[i];
        sumw += std::abs(mod_W[i]);
    }
    LOG(WARNING) << "Checking mean mode <W*dt> = " << sumw / N * mod_dt;
    LOG(WARNING) << "Checking max  mode <W*dt> = " << maxw * mod_dt;
    return 0;
}

int BO_ForceField::ForceField_write(std::ofstream& ofs0, double* nr, double* np, double* nm, const int& rdim,
                                    const int& itraj, const int& isamp) {
    // save nuclear data
    // std::cout << rdim << std::endl;
    ofs0 << FMT(8) << rdim << std::endl << FMT(8) << isamp * iou.time * mod_dt << std::endl;
    for (int j = 0; j < rdim; ++j) ofs0 << FMT(8) << nr[j] << FMT(8) << np[j] << std::endl;
    return 0;
}

Nad_ForceField::Nad_ForceField(const Param& iparm) : BO_ForceField(iparm) {
    F       = Param_GetT(int, parm, "fdim");  // Num. for electronic DOFs
    mod_occ = Param_GetT(int, parm, "occ", 0);
    // some auxiliary sizes
    FF   = F * F;
    NF   = N * F;
    NFF  = N * F * F;
    NNF  = N * N * F;
    NNFF = N * N * F * F;

    ALLOCATE_PTR_TO_VECTOR(mod_eac, F);
    ALLOCATE_PTR_TO_VECTOR(mod_rho, FF);

    if (mod_occ == -1) {
        std::string rho_read_file = Param_GetT(std::string, parm, "E_read_file");
        std::ifstream ifs(rho_read_file);
        num_real tmp;
        for (int i = 0; i < F; ++i)
            if (ifs >> tmp) mod_eac[i] = tmp;
        for (int i = 0; i < F; ++i)
            if (ifs >> tmp) mod_eac[i] += phys::math::im * tmp;
        ifs.close();
    } else if (mod_occ == -2) {
        std::string rho_read_file = Param_GetT(std::string, parm, "E_read_file");
        std::ifstream ifs(rho_read_file);
        num_real tmp;
        for (int i = 0; i < FF; ++i)
            if (ifs >> tmp) mod_rho[i] = tmp;
        for (int i = 0; i < FF; ++i)
            if (ifs >> tmp) mod_rho[i] += phys::math::im * tmp;
        ifs.close();
    }
}
Nad_ForceField::~Nad_ForceField(){};

int Nad_ForceField::get_F() { return F; }

int Nad_ForceField::ForceField_spec(double* nr, double* np, double* nm, const int& rdim, const int& fdim) {
    return (rdim > 0 && fdim > 0) ? 0 : 1;  // default, only 1 specification
}

int Nad_ForceField::nspec() { return ForceField_spec(nullptr, nullptr, nullptr, -1, -1); }

int Nad_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                    const int& rdim) {
    V[0] = 0.0f;
    for (int i = 0; i < rdim; ++i) dV[i] = 0.0f;
    return 0;
}
int Nad_ForceField::ForceField_npes(double* V, double* dV, double* ddV, double* R, double* P, const int& flag,
                                    const int& rdim, const int& itraj, const int& isamp) {
    return ForceField_npes(V, dV, ddV, R, P, flag,
                           rdim);  ///< @warning don't implement this function when you use first_call scheme
}

int Nad_ForceField::ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac, int& eocc,
                                    const int& rdim, const int& fdim, const int& itraj) {
    return 0;
}
int Nad_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                    const int& fdim) {
    return 0;
}
int Nad_ForceField::ForceField_epes(double* V, double* dV, double* ddV, double* R, const int& flag, const int& rdim,
                                    const int& fdim, const int& itraj, const int& isamp) {
    return ForceField_epes(V, dV, ddV, R, flag, rdim, fdim);
}

// Harmonic Wigner/Wavepacket sampling.
int Nad_ForceField::ForceField_init_default(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac,
                                            int& eocc, const int& rdim, const int& fdim, const int& itraj) {
    // only ask occ
    if (nr == nullptr || np == nullptr || nm == nullptr || erho == nullptr || eeac == nullptr || rdim <= 0 ||
        fdim <= 0) {
        eocc = mod_occ;
        return 0;
    }
    BO_ForceField::ForceField_init_default(nr, np, nm, rdim, itraj);
    ForceField_init_elec(erho, eeac, eocc, fdim, itraj);
    return 0;
}

int Nad_ForceField::ForceField_init_elec(num_complex* erho, num_complex* eeac, int& eocc, const int& fdim,
                                         const int& itraj) {
    eocc = mod_occ;
    if (mod_occ == -1) {
        for (int i = 0; i < F; ++i) eeac[i] = mod_eac[i];
    } else if (mod_occ == -2) {
        for (int i = 0; i < FF; ++i) erho[i] = mod_rho[i];
    }
    return 0;
}

int Nad_ForceField::CheckForceField() {
    BO_ForceField::CheckForceField();
    return 0;
}

int Nad_ForceField::ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, double* nr, double* np, double* nm,
                                     num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                     const int& itraj, const int& isamp) {
    BO_ForceField::ForceField_write(ofs0, nr, np, nm, rdim, itraj, isamp);
    // save electronic data
    ofs1 << FMT(8) << F << std::endl << FMT(8) << isamp << FMT(8) << eocc << std::endl;
    for (int i = 0, idx = 0; i < F; ++i) {
        ofs1 << FMT(8) << REAL_OF(eeac[i]) << FMT(8) << IMAG_OF(eeac[i]);
        for (int k = 0; k < F; ++k, ++idx) { ofs1 << FMT(8) << REAL_OF(erho[idx]) << FMT(8) << IMAG_OF(erho[idx]); }
        ofs1 << std::endl;
    }
    return 0;
}

int Nad_ForceField::reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim) {
    int Fadd1      = F + 1;
    num_real* pdHj = dH;
    for (int j = 0; j < N; ++j) {
        fx[j] = 0.0f;
        for (int i = 0, idxd = 0; i < F; ++i, idxd += Fadd1) {
            fx[j] += REAL_OF(rho[idxd]) * pdHj[idxd];
            for (int k = i + 1, idx1 = idxd + 1, idx2 = idxd + F; k < F; ++k, ++idx1, idx2 += F) {
                fx[j] += REAL_OF(rho[idx1]) * pdHj[idx2] + REAL_OF(rho[idx2]) * pdHj[idx1];
            }
        }
        pdHj += FF;
    }
    return 0;
}

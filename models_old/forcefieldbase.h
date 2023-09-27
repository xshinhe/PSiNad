#ifndef MODEL_BASIC_H
#define MODEL_BASIC_H

#include "model.h"

enum ForceFieldType {
    ForceFieldModel,
    ForceFieldOnTheFly,
};

class ForceField : public Model {
   public:
    ForceField(const Param& iparm);
    ForceField(const std::string& iparm_str) : ForceField(Param::parse(iparm_str)){};
    virtual ~ForceField();
    int type = ForceFieldModel;
};

class BO_ForceField : public ForceField {
   public:
    BO_ForceField(const Param& iparm);
    BO_ForceField(const std::string& iparm_str) : BO_ForceField(Param::parse(iparm_str)){};
    virtual ~BO_ForceField();
    int get_N();
    int get_Ndim();
    double Suggest_dt();
    double Suggest_tend();

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int nspec();

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp);


    int ForceField_init_default_build(const double& beta, const int& rdim);

    // Harmonic Wigner/Wavepacket sampling.
    int ForceField_init_default(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& itraj);


    virtual int CheckForceField();

    virtual int ForceField_write(std::ofstream& ofs0, num_real* nr, num_real* np, num_real* nm, const int& rdim,
                                 const int& itraj, const int& isamp);

   public:
    DEFINE_POINTER(num_real, mod_M);
    DEFINE_POINTER(num_real, mod_W);
    DEFINE_POINTER(num_real, mod_R0);
    DEFINE_POINTER(num_real, mod_P0);
    DEFINE_POINTER(num_real, mod_sigmaR);
    DEFINE_POINTER(num_real, mod_sigmaP);

   protected:
    int N, NN, Ndim;                        // for size
    num_real mod_dt, mod_tbegin, mod_tend;  // for numerical
};

class Nad_ForceField : public BO_ForceField {
   public:
    Nad_ForceField(const Param& iparm);
    Nad_ForceField(const std::string& iparm_str) : Nad_ForceField(Param::parse(iparm_str)){};
    virtual ~Nad_ForceField();
    int get_F();

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim);

    virtual int nspec();

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim, const int& itraj, const int& isamp);

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp);


    // Harmonic Wigner/Wavepacket sampling.
    int ForceField_init_default(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj);

    int ForceField_init_elec(num_complex* erho, num_complex* eeac, int& eocc, const int& fdim, const int& itraj);

    virtual int CheckForceField();

    virtual int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& isamp);

    virtual int reduce_force(num_real* fx, num_complex* rho, num_real* dH, const int& rdim, const int& fdim);

   public:
    DEFINE_POINTER(num_complex, mod_eac);
    DEFINE_POINTER(num_complex, mod_rho);

   protected:
    int F, FF, NF, NFF, NNF, NNFF;  // for size
    int mod_occ;
};

#endif  // MODEL_BASIC_H

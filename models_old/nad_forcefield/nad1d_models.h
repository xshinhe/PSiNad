#ifndef NAD1D_MODELS_H
#define NAD1D_MODELS_H

#include "../forcefieldbase.h"

namespace NAD1DPolicy {
enum _enum {
    MORSE3A,  // 3-state MORSE3A model
    MORSE3B,  // 3-state MORSE3B model
    MORSE3C,  // 3-state MORSE3C model
    MORSE15,  // 15-state MORSE model
    IVP1,     // iverted potential model 1
    IVP2,     // iverted potential model 2
    IVP3,     // iverted potential model 3
    IVP4,     // iverted potential model 4
    CL1D,     // 1-mode Caldeira-Leggett model
    JC1D,     // 1-mode Jaynes-Cummings model
    NA_I      // Na + I collision model
};
const std::map<std::string, _enum> _dict = {
    {"morse3a", MORSE3A}, {"morse3b", MORSE3B}, {"morse3c", MORSE3C}, {"morse15", MORSE15},
    {"ivp1", IVP1},       {"ivp2", IVP2},       {"ivp3", IVP3},       {"ivp4", IVP4},
    {"cl1d", CL1D},       {"jc1d", JC1D},       {"na_i", NA_I},
};
};  // namespace NAD1DPolicy

class NAD1D_ForceField : public Nad_ForceField {
   public:
    NAD1D_ForceField(const Param& iparm, const int& child);
    NAD1D_ForceField(const Param& iparm);
    NAD1D_ForceField(const std::string& iparm_str) : NAD1D_ForceField(Param::parse(iparm_str)){};


    virtual ~NAD1D_ForceField(){};

    static inline std::string name() { return "nad1d"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    virtual int NAD1D_plot();

    virtual int ForceField_epes_Morse3A(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_Morse3B(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_Morse3C(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_Morse15(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                        const int& rdim, const int& fdim);

    virtual int ForceField_epes_IVP1(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);


    virtual int ForceField_epes_IVP2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);


    virtual int ForceField_epes_IVP3(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);

    virtual int ForceField_epes_IVP4(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);


    virtual int ForceField_epes_CL1D(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);

    virtual int ForceField_epes_JC1D(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);

    virtual int ForceField_epes_NA_I(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);


   protected:
    double pm[10];

    int fftype;
    std::string ffflag;
    int spec;
};


#endif  // NAD1D_MODELS_H
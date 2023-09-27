#ifndef SPECTRUM_NAD_WRAPPER_H
#define SPECTRUM_NAD_WRAPPER_H

#include "../forcefieldbase.h"

class Spectrum_NAD_ForceField : public Nad_ForceField {
   public:
    Spectrum_NAD_ForceField(const Param& iparm);
    Spectrum_NAD_ForceField(const std::string& iparm_str) : Spectrum_NAD_ForceField(Param::parse(iparm_str)){};

    virtual ~Spectrum_NAD_ForceField();

    static inline std::string name() { return "spectrum"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    Nad_ForceField* pFF_inner;

    DEFINE_POINTER(num_real, workr_v);
    DEFINE_POINTER(num_real, workr_dv);

    bool first_call = true;
    int Fminus1;
    double ground_shift;
};


#endif  // SPECTRUM_WRAPPER_H

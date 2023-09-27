#ifndef MANYSITE_MODELS_H
#define MANYSITE_MODELS_H

#include "../forcefieldbase.h"

namespace Pauli {
const num_complex X[4] = {phys::math::iz, phys::math::iu, phys::math::iu, phys::math::iz};
const num_complex Y[4] = {phys::math::iz, -phys::math::im, phys::math::im, phys::math::iz};
const num_complex Z[4] = {phys::math::iu, phys::math::iz, phys::math::iz, -phys::math::iu};
};  // namespace Pauli

class ManySite_ForceField : public Nad_ForceField {
   public:
    // enum HamType { Landau_Zener, DEMKOV_OSHEROV, BOW_TIE, Anderson, Hubbard, Ising, XYModel };
    // ManySite_ForceField(const Param& iparm, const int& child);

    ManySite_ForceField(const Param& iparm);
    ManySite_ForceField(const std::string& iparm_str) : ManySite_ForceField(Param::parse(iparm_str)){};


    virtual ~ManySite_ForceField(){};

    static inline std::string name() { return "manysite"; }

    virtual int ForceField_heff(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim);

    virtual int ForceField_heff_Ising(num_complex* H, num_complex* rhos, const int& mdim, const int& fdim);

    inline int get_M() { return M; }

    int M;  // M sites

    DEFINE_POINTER(num_real, JpMat);
    DEFINE_POINTER(num_real, JzMat);
    DEFINE_POINTER(num_complex, redX);
    DEFINE_POINTER(num_complex, redY);
    DEFINE_POINTER(num_complex, redZ);

    double Jp, Jz, alpha, omega;
};

#endif  // MANYSITE_MODELS_H
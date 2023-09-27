#ifndef ZnPc_H
#define ZnPc_H

#include "systembath.h"

class ZnPc_ForceField : public SystemBath_ForceField {
   public:
    ZnPc_ForceField(const Param& iparm);
    ZnPc_ForceField(const std::string& iparm_str) : ZnPc_ForceField(Param::parse(iparm_str)){};
    virtual ~ZnPc_ForceField(){};

    static inline std::string name() { return "znpc"; }

    int index(const int& iL, const int& ie, const int& ih);

    int init_Hamiltonian();

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim, const int& itraj, const int& isamp);

    virtual int get_nbath();
    virtual int get_Nb();

   protected:
    int nmol, nexc;

    DEFINE_POINTER_PROTECTED(int, idxarr_L);
    DEFINE_POINTER_PROTECTED(int, idxarr_es);
    DEFINE_POINTER_PROTECTED(int, idxarr_hs);

    DEFINE_POINTER_PROTECTED(num_real, dQe1);
    DEFINE_POINTER_PROTECTED(num_real, dQe2);
    DEFINE_POINTER_PROTECTED(num_real, dQc);
    DEFINE_POINTER_PROTECTED(num_real, dQa);
    DEFINE_POINTER_PROTECTED(num_real, w2dQe1);
    DEFINE_POINTER_PROTECTED(num_real, w2dQe2);
    DEFINE_POINTER_PROTECTED(num_real, w2dQc);
    DEFINE_POINTER_PROTECTED(num_real, w2dQa);

    DEFINE_POINTER_PROTECTED(num_real, Etilde);
    DEFINE_POINTER_PROTECTED(num_real, Vtilde);
    DEFINE_POINTER_PROTECTED(num_real, te_tilde);
    DEFINE_POINTER_PROTECTED(num_real, th_tilde);
    DEFINE_POINTER_PROTECTED(num_real, tect_tilde);
    DEFINE_POINTER_PROTECTED(num_real, thct_tilde);
    DEFINE_POINTER_PROTECTED(num_real, eigen_E);
    DEFINE_POINTER_PROTECTED(num_real, eigen_T);
};

#endif  // ZnPc_H

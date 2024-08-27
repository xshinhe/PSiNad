#ifndef Model_ZnPc_H
#define Model_ZnPc_H

#include "systembath.h"

class Model_ZnPc : public Model {
   public:
   protected:
    int             nmol, nexc;
    span<kids_int>  idxarr_L;
    span<kids_int>  idxarr_es;
    span<kids_int>  idxarr_hs;
    span<kids_real> dQe1;
    span<kids_real> dQe2;
    span<kids_real> dQc;
    span<kids_real> dQa;
    span<kids_real> w2dQe1;
    span<kids_real> w2dQe2;
    span<kids_real> w2dQc;
    span<kids_real> w2dQa;
    span<kids_real> Etilde;
    span<kids_real> Vtilde;
    span<kids_real> te_tilde;
    span<kids_real> th_tilde;
    span<kids_real> tect_tilde;
    span<kids_real> thct_tilde;
    span<kids_real> eigen_E;
    span<kids_real> eigen_T;
};

#endif  // Model_ZnPc_H

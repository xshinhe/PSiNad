#ifndef Model_ZnPc_H
#define Model_ZnPc_H

#include "systembath.h"

class Model_ZnPc : public Model {
   public:
   protected:
    int             nmol, nexc;
    span<psnd_int>  idxarr_L;
    span<psnd_int>  idxarr_es;
    span<psnd_int>  idxarr_hs;
    span<psnd_real> dQe1;
    span<psnd_real> dQe2;
    span<psnd_real> dQc;
    span<psnd_real> dQa;
    span<psnd_real> w2dQe1;
    span<psnd_real> w2dQe2;
    span<psnd_real> w2dQc;
    span<psnd_real> w2dQa;
    span<psnd_real> Etilde;
    span<psnd_real> Vtilde;
    span<psnd_real> te_tilde;
    span<psnd_real> th_tilde;
    span<psnd_real> tect_tilde;
    span<psnd_real> thct_tilde;
    span<psnd_real> eigen_E;
    span<psnd_real> eigen_T;
};

#endif  // Model_ZnPc_H

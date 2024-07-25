#ifndef Model_ZnPc_H
#define Model_ZnPc_H

#include "systembath.h"

class Model_ZnPc : public Model {
   public:
   protected:
    int        nmol, nexc;
    kids_int*  idxarr_L;
    kids_int*  idxarr_es;
    kids_int*  idxarr_hs;
    kids_real* dQe1;
    kids_real* dQe2;
    kids_real* dQc;
    kids_real* dQa;
    kids_real* w2dQe1;
    kids_real* w2dQe2;
    kids_real* w2dQc;
    kids_real* w2dQa;
    kids_real* Etilde;
    kids_real* Vtilde;
    kids_real* te_tilde;
    kids_real* th_tilde;
    kids_real* tect_tilde;
    kids_real* thct_tilde;
    kids_real* eigen_E;
    kids_real* eigen_T;
};

#endif  // Model_ZnPc_H

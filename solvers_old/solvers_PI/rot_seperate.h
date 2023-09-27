#ifndef ROT_SEPERATE_H
#define ROT_SEPERATE_H

#include "../../utils/array_eigen.h"
#include "../../utils/types.h"
#include
void rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in);
void pseudo_inv(int N, num_real* A, num_real* invA, num_real eps);
void cross(num_real* vec1, num_real* vec2, num_real* prod);
num_real *idx, *idp, *idf, *idm;  // temporary cursor variables
num_real *xr, *vr, *fr;           // centroid relative variables
num_real Xc[3] = {0, 0, 0};
num_real Vc[3] = {0, 0, 0};
num_real Fc[3] = {0, 0, 0};
num_real Lc[3] = {0, 0, 0};
num_real Wc[3] = {0, 0, 0};
num_real Mc[3] = {0, 0, 0};
num_real Ac[3] = {0, 0, 0};
num_real vectmp[3];
num_real Ic[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
num_real invIc[9];
num_real Mt;
#endif
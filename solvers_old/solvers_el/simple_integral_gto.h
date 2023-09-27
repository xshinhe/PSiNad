#ifndef Simple_Integral_GTO_H
#define Simple_Integral_GTO_H

#include <cmath>  // erf
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "../../utils/phys.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
/*

        the general integral function as follow
        all integral formula refer <<Modern quantum chemistry computational methods>>, Zhizhong Wang.
        and << quantum chenistry >> by Guangxian Xu, Lemin Li.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////


int comb(const int& a, const int& b);
int fact(const int& n);
int dfact(const int& n);

inline double gto_normfactor(const double& alpha, const int& l, const int& m, const int& n) {
    int lmn = l + m + n;
    return sqrt(pow(2, 2 * lmn) * pow(alpha, lmn + 1.5) * pow(2 / phys::math::pi, 1.5) /
                (dfact(2 * l - 1) * dfact(2 * m - 1) * dfact(2 * n - 1)));
}

inline int parity(const int& n) { return (n % 2 == 0 ? 1 : -1); }


double fun_f(const int& l, const int& l1, const int& l2, const double& a, const double& b);

double fun_G(const int& I, const int& l1, const int& l2, const double& a, const double& b, const double& c,
             const double& w);


double fun_H(const int& l, const int& l1, const int& l2, const double& a, const double& b, const double& w);

double fun_D(const int& I, const int& l1, const int& l2, const int& l3, const int& l4, const double& PQ,
             const double& a, const double& b, const double& c, const double& d, const double& del, const double& a1,
             const double& a2, const double& a3, const double& a4);
double gto_S(const double& a1, const double& x1, const double& y1, const double& z1, const int& l1, const int& m1,
             const int& n1, const double& a2, const double& x2, const double& y2, const double& z2, const int& l2,
             const int& m2, const int& n2);

double gto_T(const double& a1, const double& x1, const double& y1, const double& z1, const int& l1, const int& m1,
             const int& n1, const double& a2, const double& x2, const double& y2, const double& z2, const int& l2,
             const int& m2, const int& n2);

double gto_V(const double& a1, const double& x1, const double& y1, const double& z1, const int& l1, const int& m1,
             const int& n1, const double& a2, const double& x2, const double& y2, const double& z2, const int& l2,
             const int& m2, const int& n2, const double& xc, const double& yc, const double& zc);

double gto_ERI(const double& a1, const double& x1, const double& y1, const double& z1, const int& l1, const int& m1,
               const int& n1, const double& a2, const double& x2, const double& y2, const double& z2, const int& l2,
               const int& m2, const int& n2, const double& a3, const double& x3, const double& y3, const double& z3,
               const int& l3, const int& m3, const int& n3, const double& a4, const double& x4, const double& y4,
               const double& z4, const int& l4, const int& m4, const int& n4);

double Fn(const int& m, const double& w);


//////////////////////////////////////////////////////////////////////////////////////////////////
// declarition of functions

double IntecGTO_S(double* clist1, double* alist1, const double& x1, const double& y1, const double& z1, const int& L1,
                  const int& M1, const int& N1, const int& nc1, double* clist2, double* alist2, const double& x2,
                  const double& y2, const double& z2, const int& L2, const int& M2, const int& N2, const int& nc2);

double IntecGTO_T(double* clist1, double* alist1, const double& x1, const double& y1, const double& z1, const int& L1,
                  const int& M1, const int& N1, const int& nc1, double* clist2, double* alist2, const double& x2,
                  const double& y2, const double& z2, const int& L2, const int& M2, const int& N2, const int& nc2);

double IntecGTO_V(double* clist1, double* alist1, const double& x1, const double& y1, const double& z1, const int& L1,
                  const int& M1, const int& N1, const int& nc1, double* clist2, double* alist2, const double& x2,
                  const double& y2, const double& z2, const int& L2, const int& M2, const int& N2, const int& nc2,
                  const double& xc, const double& yc, const double& zc, const int& znum);

double IntecGTO_ERI(double* clist1, double* alist1, const double& x1, const double& y1, const double& z1, const int& L1,
                    const int& M1, const int& N1, const int& nc1, double* clist2, double* alist2, const double& x2,
                    const double& y2, const double& z2, const int& L2, const int& M2, const int& N2, const int& nc2,
                    double* clist3, double* alist3, const double& x3, const double& y3, const double& z3, const int& L3,
                    const int& M3, const int& N3, const int& nc3, double* clist4, double* alist4, const double& x4,
                    const double& y4, const double& z4, const int& L4, const int& M4, const int& N4, const int& nc4);

//////////////////////////////////////////////////////////////////////////////////////////////////


#endif  // Simple_Integral_GTO_H

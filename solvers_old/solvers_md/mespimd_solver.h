#ifndef SOLVER_MESPIMD_H
#define SOLVER_MESPIMD_H

#include "pimd_solver.h"

namespace MESPIMDFactorizePolicy {
enum _enum { Diag, Hyper, First };
const std::map<std::string, _enum> _dict = {{"diag", Diag}, {"hyper", Hyper}, {"first", First}};
};  // namespace MESPIMDFactorizePolicy

class MESPIMDTraj_Solver : public PIMDTraj_Solver {
   public:
    enum mespimd_enum { firt_order, hyperbolic, diagonalization };

    MESPIMDTraj_Solver(const Param& iparm, Model* pM);
    MESPIMDTraj_Solver(const std::string& iparm_str, Model* pM) : MESPIMDTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~MESPIMDTraj_Solver();

    static inline std::string name() { return "mespimd"; }

    virtual int init(const int& itraj);

    virtual int ff_calc1(const int& level = 1);

    virtual int ff_Oi(int i, num_real* Oi_pos, num_real* dVi_pos, num_real* Oi, num_real* dOi, num_real* Vi,
                      num_real* dVi);

    virtual int ff_OO();

    virtual num_real esti_term1(num_real* Q, const bool& dressed = true, const bool& fixed = false);

    virtual num_real esti_term2(num_real* Q1, num_real* Q2, const bool& dressed1 = true, const bool& dressed2 = true,
                                const bool& fixed1 = false, const bool& fixed2 = false);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

   protected:
    int F, FF, NF, NFF, NNF, NNFF, FFFF;

    DEFINE_POINTER_PROTECTED(num_real, Vs);
    DEFINE_POINTER_PROTECTED(num_real, dVs);
    DEFINE_POINTER_PROTECTED(num_real, ddVs);
    DEFINE_POINTER_PROTECTED(num_real, Es);
    DEFINE_POINTER_PROTECTED(num_real, dEs);
    DEFINE_POINTER_PROTECTED(num_real, Ts);
    DEFINE_POINTER_PROTECTED(num_real, mat_fA);  // exact potential part of partition function matrix
    DEFINE_POINTER_PROTECTED(num_real, mat_fD);  // positive-defined potential part of partition function matrix

    num_real tr_fA, tr_fD, rw_factor;

    DEFINE_POINTER_PROTECTED(num_real, O_pos);
    DEFINE_POINTER_PROTECTED(num_real, dV_pos);
    DEFINE_POINTER_PROTECTED(num_real, O);
    DEFINE_POINTER_PROTECTED(num_real, OO);      ///< OO = O^* O
    DEFINE_POINTER_PROTECTED(num_real, OOb);     ///< \prod_{k<b} OO(k)
    DEFINE_POINTER_PROTECTED(num_real, OOe);     ///< \prod_{k>e} OO(k)
    DEFINE_POINTER_PROTECTED(num_real, OObe);    ///< \prod_{b<k<e} OO(k)
    DEFINE_POINTER_PROTECTED(num_real, dO);      ///< dO = \partial O / \patial_R
    DEFINE_POINTER_PROTECTED(num_real, V2);      ///< V2 = V * V
    DEFINE_POINTER_PROTECTED(num_real, hRdOO);   ///< note: hRd(...) = 0.5 * (R) * \partial_R(...)
    DEFINE_POINTER_PROTECTED(num_real, hRcdOO);  ///< note: hRcd(...) = 0.5 * (R-Rc) * \partial_R(...)
    DEFINE_POINTER_PROTECTED(num_real, hRdOVO);
    DEFINE_POINTER_PROTECTED(num_real, hRcdOVO);
    DEFINE_POINTER_PROTECTED(num_real, rho);
    DEFINE_POINTER_PROTECTED(num_real, rho_op);

    // only used for calling ForceField_init
    int occ0;
    DEFINE_POINTER_PROTECTED(num_complex, eac0);
    DEFINE_POINTER_PROTECTED(num_complex, rho0);

    Nad_ForceField* pForceField;

    std::string factor_flag;
    int factor_type;
};

#endif  // SOLVER_MESPIMD_H

#ifndef PIMDTraj_Solver_H
#define PIMDTraj_Solver_H

#include "../../utils/la_utils.h"
#include "traj.h"

namespace PIMDTransformPolicy {
enum _enum { Primitive, Staging, NormalMode, TEST };
const std::map<std::string, _enum> _dict = {
    {"#prim", Primitive}, {"#stag", Staging}, {"#norm", NormalMode}, {"#test", TEST}};
};  // namespace PIMDTransformPolicy

namespace PIMDIntegratorPolicy {
enum _enum { BAOAB, BCOCB };
const std::map<std::string, _enum> _dict = {{"BAOAB", BAOAB}, {"BCOCB", BCOCB}};
};  // namespace PIMDIntegratorPolicy

namespace PIMDRestartPolicy {
enum _enum { none, in, out, both };
const std::map<std::string, _enum> _dict = {{"none",none}, {"in",in}, {"out",out}, {"both",both}};
};  
/**
 * PIMDTraj_solver is enhanced from Traj_Solver
 */
class PIMDTraj_Solver : public Traj_Solver {
   public:
    PIMDTraj_Solver(Param iparm, Model* pM);
    PIMDTraj_Solver(const std::string& iparm_str, Model* pM) : PIMDTraj_Solver(Param::parse(iparm_str), pM){};

    virtual ~PIMDTraj_Solver();

    static inline std::string name() { return "pimd"; }

    virtual int ff_calc1(const int& level);

    virtual int spring_force();

    virtual int update_r(const num_real& dt);

    virtual int update_p(const num_real& dt);

    virtual int update_p_harm(const num_real& dt);

    virtual int caylay_update_half(
        const num_real& dt);  // half step of cayley modification is cay(AΔt)^0.5 rather than cay(0.5AΔt)

    virtual int BAOAB(int& succ,
                      int step);  // BAOAB-num in J. Chem. Phys. 145, 024103 (2016) https://doi.org/10.1063/1.4954990

    virtual int BCOCB(int& succ, int step);

    virtual int update_thermo(const num_real& dt);

    virtual int traj(TCFnucl& tcfer, const int& PN);

    virtual int traj_Middle(TCFnucl& tcfer, const int& PN);

    virtual int traj_property(const num_real& dt);

    virtual int sampler(const int& isamp, TCFnucl& tcfer);

    virtual int estimator(const int& isamp, TCFnucl& tcfer);

    virtual int run_impl();

    virtual int run_parallel();

    virtual int init(const int& itraj);

    int all_X2K();

    int all_K2X();

    int all_FX2FK();

    virtual int rot_trans_corr(int Natom, num_real* m_in, num_real* x_in, num_real* p_in, num_real* F_in,
                               bool cal_force = true);

    virtual int pseudo_inv(int N, num_real* A, num_real* invA, num_real* vectmp, num_real eps);

    virtual int cross(num_real* vec1, num_real* vec2, num_real* prod);

    virtual int rst_output(const int& traj_in);

    virtual int rst_read(const int& traj_in);

   protected:
    int P, PP, PN, PNPN, Natom;
    int pimd_type, integ_type, rst_type;
    num_real bf2, betap;

    DEFINE_POINTER_PROTECTED(num_real, nrs);   ///< real position for PIMD (old nr is temporarily used ad centroid bead)
    DEFINE_POINTER_PROTECTED(num_real, nks);   ///< transformed position for PIMD
    DEFINE_POINTER_PROTECTED(num_real, nps);   ///< virtual momentum for PIMD (conjugated with nks, not nrs)
    DEFINE_POINTER_PROTECTED(num_real, nms);   ///< virtual mass
    DEFINE_POINTER_PROTECTED(num_real, nfs);   ///< forces for PIMD
    DEFINE_POINTER_PROTECTED(num_real, nfks);  ///< transformed forces for PIMD
    DEFINE_POINTER_PROTECTED(num_real, masswgt);  ///< arbitrary weighting factors for mass
    DEFINE_POINTER_PROTECTED(num_real, bfwgt);  ///< modified weighting factors for frequency (based on transformation)
    DEFINE_POINTER_PROTECTED(num_real, D2);     ///< spring interaction matrix
    DEFINE_POINTER_PROTECTED(num_real, Tran);   ///< transform matrix between nrs and nks
    DEFINE_POINTER_PROTECTED(num_real, fk_spring);   ///< spring force of pimd of transformed coordinates

    DEFINE_POINTER_PROTECTED(num_real, vpeses);
    DEFINE_POINTER_PROTECTED(num_real, grads);
    DEFINE_POINTER_PROTECTED(num_real, hesses);

    num_real *rc, *vc, *ac, *Lc, *Wc, *Mc, *Ac, *vectmp, *Mattmp, *Ic, *invIc, Mt;
};


#endif  // PIMDTraj_Solver_H

/**
 * @file Kernel_Elec.h
 * @author xshinhe
 * @version 1.1
 * @date 2023-03
 * @brief initialization kernels for electonic DOFs
 * @details
 *  The initialization of electonic DOFs are tightly related to Solver.
 *  Use it in Solver's Kernel_Builder();
 */

#ifndef Kernel_Elec_CMSH_H
#define Kernel_Elec_CMSH_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

DEFINE_POLICY(CMSHPolicy,
              EHR,    // Ehrenfest Dynamics
              BOSH,   // BO dynamics & hopping according to W(\rho)
              CVSH,   // CV dynamics & hopping acoording to W(\rho)
              BOSD,   // BO dynamics & smoothing acoording to W(\rho)
              CVSD);  // CV dynamics & smoothing acoording to W(\rho)

/**
 * @brief initialization kernel for electonic DOFs in CMSH
 */
class Kernel_Elec_CMSH final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_CMSH"; }

    Kernel_Elec_CMSH(double scale = 1.0e0) : Kernel(), scale{scale} {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

   private:
    CMSHPolicy::_type cmsh_type;

    kids_real gamma1, gamma2, xi1, xi2;
    bool use_focus = false;
    bool use_cv    = true;   // adapt cv in rho_nuc
    bool use_wmm   = false;  // in this case, gamma1 will be used as delta in wMM
    bool use_fall  = false;
    bool use_gdtwa = false;
    bool use_sum   = false;

    bool cread_from_ds        = false;
    bool disable_inner_switch = false;

    bool reflect = true;  // treatment in hopping
    int hopping_type1;
    int hopping_type2;

    bool dynamic_alpha;

    double scale;
    double dt;
    kids_real alpha0;
    kids_real* alpha;
    kids_real* p;
    kids_real* m;
    kids_real* fadd;
    kids_real* ftmp;
    kids_real* direction;
    kids_real *vpes, *V, *E, *dE, *T;
    kids_real* Epot;
    kids_complex* H;
    kids_complex* wrho;
    kids_real* sqcw;
    kids_real *sqcIA, *sqcID;

    bool* do_prec_ptr;

    //
    bool use_sqc;
    int sqc_init;
    bool only_adjust;
    bool check_cxs;

    bool use_fssh;
    bool use_strange_win;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_CMSH_H

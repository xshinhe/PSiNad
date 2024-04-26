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

#ifndef Kernel_Elec_Switch_H
#define Kernel_Elec_Switch_H

#include "kids/Kernel.h"
#include "kids/Kernel_Elec.h"
#include "kids/Kernel_Elec_NAD.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in NAD
 */
class Kernel_Elec_Switch final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Kernel_Elec_Switch(double scale = 1.0e0) : Kernel(), scale{scale} {
        appendChild(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

   private:
    NADPolicy::_type cmsh_type;

    kids_real gamma1, gamma2, xi1, xi2;
    bool      use_focus = false;
    bool      use_cv    = true;   // adapt cv in rho_nuc
    bool      use_wmm   = false;  // in this case, gamma1 will be used as delta in wMM
    bool      use_fall  = false;
    bool      use_gdtwa = false;
    bool      use_sum   = false;

    bool cread_from_ds = false;

    bool reflect = true;  // treatment in hopping
    int  hopping_type1;
    int  hopping_type2;

    bool dynamic_alpha;

    double        scale;
    double        dt;
    kids_real     alpha0;
    kids_real*    alpha;
    kids_real*    p;
    kids_real*    m;
    kids_real*    fadd;
    kids_real*    ftmp;
    kids_real*    direction;
    kids_real *   vpes, *V, *E, *dE, *T;
    kids_real*    Epot;
    kids_complex* H;
    kids_complex* wrho;
    kids_real*    sqcw;
    kids_real *   sqcIA, *sqcID;

    kids_bint* at_samplingstep_finally_ptr;

    //
    bool use_sqc;
    int  sqc_init;
    bool only_adjust;
    bool check_cxs;

    bool use_fssh;
    bool use_strange_win;

    virtual void setInputParam_impl(std::shared_ptr<Param>& PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_Switch_H

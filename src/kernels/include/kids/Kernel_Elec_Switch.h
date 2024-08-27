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

namespace PROJECT_NS {

/**
 * @brief initialization kernel for electonic DOFs in NAD
 */
class Kernel_Elec_Switch final : public Kernel {
   public:
    virtual const std::string getName();

    virtual int getType() const;

   private:
    bool reflect = true;  // treatment in hopping
    int  hopping_direction_type;
    int  hopping_choose_type;

    span<kids_int>     occ_nuc;
    span<kids_real>    dt_ptr;
    span<kids_real>    T;
    span<kids_real>    Epot;
    span<kids_real>    vpes;
    span<kids_real>    p;
    span<kids_real>    m;
    span<kids_real>    EMat;
    span<kids_real>    ForceMat;
    span<kids_real>    direction;
    span<kids_complex> rho_ele;
    span<kids_complex> rho_nuc;
    span<kids_complex> H;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Elec_Switch_H

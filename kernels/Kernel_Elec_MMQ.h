#ifndef Kernel_Elec_MMQ_H
#define Kernel_Elec_MMQ_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

DEFINE_POLICY(MMQPolicy,  //
              DUAL);

/**
 * @brief initialization kernel for electonic DOFs in MMQ
 */
class Kernel_Elec_MMQ final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_MMQ"; }

    Kernel_Elec_MMQ() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

   private:
    MMQPolicy::_type mmq_type;
    num_real scale;
    num_real xi, gamma;
    int Fref;
    bool pure_phase;
    bool cont_phase;
    bool rand_act;

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // Kernel_Elec_MMQ_H

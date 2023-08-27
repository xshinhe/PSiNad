#ifndef Kernel_Elec_MMD_H
#define Kernel_Elec_MMD_H

#include "../core/Kernel.h"
#include "../core/Policy.h"
#include "Kernel_Elec.h"

namespace PROJECT_NS {

DEFINE_POLICY(MMDPolicy,  //
              MMF,        //
              TWA,        //
              MID);

/**
 * @brief initialization kernel for electonic DOFs in MMD
 */
class Kernel_Elec_MMD final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Elec_MMD"; }

    Kernel_Elec_MMD() {
        push(std::shared_ptr<Kernel_Elec>(new Kernel_Elec()));  //
    }

    static int rho_focus(num_complex *rho, int iocc, double gamma_ou, double gamma_uu, int fdim, bool rand_act = false,
                         bool pure_phase = true, bool cont_phase = true);

   private:
    MMDPolicy::_type mmd_type;
    num_real scale;
    num_real gamma_ou, gamma_uu;
    int Fref;
    bool pure_phase;
    bool cont_phase;
    bool rand_act;

    virtual void read_param_impl(Param *PM);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS

#endif  // Kernel_Elec_MMD_H

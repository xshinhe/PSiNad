#ifndef Kernel_Conserve_H
#define Kernel_Conserve_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * This class implements a process for energy conservation.
 */
class Kernel_Conserve final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Conserve"; }

   private:
    kids_real* E;          ///< adiabatic energies
    kids_real* p;          ///< nuclear momemtum
    kids_real* m;          ///< nuclear mass
    kids_real* Etot;       ///< total energy
    kids_real* Etot_prev;  ///< total energy in previous step
    kids_real* Etot_init;  ///< total energy at initial time
    kids_real* Ekin;       ///< kinematic energy
    kids_real* Epot;       ///< potential energy
    kids_real* vpes;       ///< potential energy for only nuclear part

    int conserve_direction;

    bool conserve_scale;
    bool* succ_ptr;
    bool* frez_ptr;
    bool* last_attempt_ptr;
    int* fail_type_ptr;

    int cnt_loose = 0;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Conserve_H

#ifndef INTERF_SANDER_H
#define INTERF_SANDER_H

#include <unistd.h>

#include <cmath>

#include "../../utils/elements.h"
#include "../forcefieldbase.h"
#include "sander.h"

class sander_ForceField : public BO_ForceField {
   public:
    sander_ForceField(const Param& iparm);

    virtual ~sander_ForceField();

    static inline std::string name() { return "sander"; }
    /**
     * @brief ForceField_init for mndo99
     * @param
     *     nr:  nuclear configuration
     *     np:  nuclear momentum
     *     nm:  nuclear mass
     *   rdim:  nuclear freedoms (=3*natom)
     * @return int
     * @bug none
     */
    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    /**
     * @brief override ForceField_epes within mndo99 (interface) and save temporary calculations
     * @param
     *      E:  adiabatic potential surface
     *     dE:  gradients at adiabatic surface & coupled gradients between adiabatic surfaces
     *    ddE:  hessian in adiabatic surfaces
     *      R:  input configuration, in order (x1,y1,z1,...,xN,yN,zN), in atomic units
     *   flag:  = 0 only calculate E
     *          = 1 calculate E and dE
     *          = 2 calculate E, dE and ddE (not supported)
     *   rdim:  nuclear freedoms (=3*natom)
     *   fdim:  electronic freedoms
     * @return int
     * @retval -1, failed
     * @retval >=0, parsed E
     * @retval >=1, parsed dE
     * @bug none
     */
    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

   protected:
    int natom, nmol;
    int read_flag;

    sander_input inp;

    std::string inpcrd, prmtop;
    pot_ene ene;
    num_real *crd, *frc;

    // std::string savefile_traj, savefile_ener, savefile_grad, savefile_nac;
    // std::string savefile_RP, savefile_eac;
};

#endif  // INTERF_SANDER_H

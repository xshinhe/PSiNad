#ifndef INTERF_OPENMM_H
#define INTERF_OPENMM_H

#include <unistd.h>

#include <cmath>

#include "../forcefieldbase.h"

class OPENMM_ForceField : public BO_ForceField {
   public:
    OPENMM_ForceField(const Param& iparm);
    OPENMM_ForceField(const std::string& iparm_str) : OPENMM_ForceField(Param::parse(iparm_str)){};

    virtual ~OPENMM_ForceField();

    /**
     * @brief ForceField_init for mndo99
     * @param
     *     nr:  nuclear configuration
     *     np:  nuclear momentum
     *     nm:  nuclear mass
     *   erho:  electronic density
     *   eeac:  electronic amplititude
     *   eocc:  electronic occupation
     *   rdim:  nuclear freedoms (=3*natom)
     *   fdim:  electronic freedoms
     *  itraj:  index of trajectory
     * @return int
     * @bug none
     */
    virtual int ForceField_init(double* nr, double* np, double* nm, num_complex* erho, num_complex* eeac, int& eocc,
                                const int& rdim, const int& fdim, const int& itraj);

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
     *  itraj:  #-th of trajectories
     *  istep:  #-th of steps
     * @return int
     * @retval -1, failed
     * @retval >=0, parsed E
     * @retval >=1, parsed dE
     * @bug none
     */
    virtual int ForceField_npes(num_real* E, num_real* dE, num_real* ddE,
                                num_real* R,                       // input in au
                                const int& flag,                   // flag
                                const int& rdim, const int& fdim,  // size
                                const int& itraj, const int& istep);

    /**
     * @brief override ForceField_write within mndo99
     * @param
     *   ofs0:  output stream0, output nuclear information
     *   ofs1:  output stream1, output electronic information, (if save to the same file, let ofs1 equals ofs0)
     *     nr:  nuclear configuration
     *     np:  nuclear momentum
     *     nm:  nuclear mass
     *   erho:  electronic density
     *   eeac:  electronic amplititude
     *   eocc:  electronic occupation
     *   rdim:  nuclear DOFs
     *   fdim:  electronic DOFs
     *  itraj:  index of trajectory
     *  istep:  index of steps
     * @return int
     * @bug none
     */
    virtual int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, double* nr, double* np, double* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& istep);

    int parse_mndo99(const std::string& mndo99inp);
    std::string new_keyword(const std::map<std::string, int>& newkeyword);
    int new_task(num_real* R, const std::string& file, const std::string& task_flag, const int& rdim);
    int track_nac_sign();
    int parse_standard(const std::string& log);
    int parse_hessian(const std::string& log);
    int calc_normalmode(num_real* R, const int& rdim);
    int calc_samp();
    int calc_scan();

   protected:
    int natom, *atoms;
    int read_flag;

    std::string savefile_traj, savefile_ener, savefile_grad, savefile_nac;
    std::string savefile_RP, savefile_eac;
};


#endif  // INTERF_OPENMM_H

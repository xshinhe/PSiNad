/**
 * @file interf_mndo99mrci.h
 * @author xshinhe
 * @version 1.0
 * @date 2021-12-31
 * @brief forcefield interface to mndo99
 * @details none
 * @see Nad_ForceField class
 * @note none
 * @todo
 */

#ifndef INTERF_MNDO99MRCI_H
#define INTERF_MNDO99MRCI_H

#include <unistd.h>

#include <cmath>

#include "../forcefieldbase.h"

struct MNDO99KW {
    std::string key;
    int val;
};

using MNDO99KW_map = std::map<std::string, int>;

class MNDO99_ForceField : public Nad_ForceField {
   public:
    MNDO99_ForceField(const Param& iparm);
    MNDO99_ForceField(const std::string& iparm_str) : MNDO99_ForceField(Param::parse(iparm_str)){};

    virtual ~MNDO99_ForceField(){};

    static inline std::string name() { return "mndo99"; }

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
    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj);


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
    virtual int ForceField_epes(num_real* E, num_real* dE, num_real* ddE,
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
    virtual int ForceField_write(std::ofstream& ofs0, std::ofstream& ofs1, num_real* nr, num_real* np, num_real* nm,
                                 num_complex* erho, num_complex* eeac, int& eocc, const int& rdim, const int& fdim,
                                 const int& itraj, const int& istep);

    int parse_mndo99(const std::string& mndo99inp);
    std::string new_keyword(const MNDO99KW_map& newkeyword);
    int new_task(num_real* R, const std::string& file, const std::string& task_flag, const int& rdim);
    int track_nac_sign();
    int parse_standard(const std::string& log);
    int parse_hessian(const std::string& log);
    int calc_normalmode(num_real* R, const int& rdim);
    int calc_samp();
    int calc_scan();

   protected:
    int natom;
    int read_flag;
    int ncigrd;
    int iroot;

    DEFINE_POINTER_PROTECTED(int, atoms);
    DEFINE_POINTER_PROTECTED(num_real, mod_Hess);
    DEFINE_POINTER_PROTECTED(num_real, mod_Tmat);
    DEFINE_POINTER_PROTECTED(num_real, nr_samp);
    DEFINE_POINTER_PROTECTED(num_real, np_samp);
    DEFINE_POINTER_PROTECTED(num_real, nac_prev);
    DEFINE_POINTER_PROTECTED(num_real, nac);
    DEFINE_POINTER_PROTECTED(num_real, ener);
    DEFINE_POINTER_PROTECTED(num_real, ener2);
    DEFINE_POINTER_PROTECTED(num_real, grad);

    num_real temp, beta;
    bool diff_nac;
    std::string init_nuclinp, savename;
    // std::vector<std::string> mndo99_keyword, mndo99_comment, mndo99_data, mndo99_addition;
    std::vector<std::string> mndo99_data;
    std::string mndo99_keyword, mndo99_comment, mndo99_addition;

    std::vector<MNDO99KW> keyword;  // keyword wrapper

    std::string savefile_traj, savefile_ener, savefile_grad, savefile_nac;
    std::string savefile_RP, savefile_eac;
};


#endif  // INTERF_MNDO99MRCI_H

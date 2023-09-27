#ifndef INTERF_GAUSSTDDFT_H_H
#define INTERF_GAUSSTDDFT_H_H

#include <unistd.h>

#include <algorithm>

#include "../forcefieldbase.h"

class GAUSS16_ForceField : public Nad_ForceField {
   public:
    GAUSS16_ForceField(const Param& iparm);
    GAUSS16_ForceField(const std::string& iparm_str) : GAUSS16_ForceField(Param::parse(iparm_str)){};

    virtual ~GAUSS16_ForceField(){};

    static inline std::string name() { return "gauss16"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_epes(num_real* E, num_real* dE, num_real* ddE,
                                num_real* R,  // input in au
                                const int& flag, const int& rdim, const int& fdim);
    int parse_g16(const std::string& g16inp);
    int calc_hess(num_real* R, const int& rdim);

   protected:
    int natom;
    int read_flag;

    DEFINE_POINTER_PROTECTED(int, atoms);
    DEFINE_POINTER_PROTECTED(num_real, mod_Hess);
    DEFINE_POINTER_PROTECTED(num_real, mod_Tmat);
    DEFINE_POINTER_PROTECTED(num_real, nr_samp);
    DEFINE_POINTER_PROTECTED(num_real, np_samp);
    num_real temp, beta;
    std::string init_nuclinp, savename;
    std::vector<std::string> g16_keyword, g16_comment, g16_data, g16_addition;
    Param keyword;  // jcson wrapper
};

#endif  // INTERF_GAUSSTDDFT_H_H

#ifndef INTERF_FF_PYTHON_H
#define INTERF_FF_PYTHON_H

#include "../forcefieldbase.h"

class PythonFF_ForceField : public Nad_ForceField {
   public:
    PythonFF_ForceField(const Param& iparm);

    PythonFF_ForceField(const std::string& iparm_str) : PythonFF_ForceField(Param::parse(iparm_str)){};

    virtual ~PythonFF_ForceField();

    static inline std::string name() { return "pyff"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& itraj);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

   protected:
    std::string pymod_file;
};

#endif

#ifndef PYRINCAVITY_H
#define PYRINCAVITY_H

#include "../forcefieldbase.h"

namespace PyrCavPolicy {
enum _enum {
    PC1,  // classical bath, Qc, Qt and cav-modes, but quantum level
          /**
           *
           * Notations:
           *      |g,(ql,qc,qt,q)>, |e1,(ql,qc,qt,q)>, |e2,(ql,qc,qt,q)>
           */
    PC2,  // classical bath, Qc, Qt, but quantum cav-modes & level
          /**
           *
           * Notations:
           *      |g,0 (qc,qt,q)>, |e1,0 (qc,qt,q)>, |e2,0 (qc,qt,q)>
           *      |g,1 (qc,qt,q)>, |e1,1 (qc,qt,q)>, |e2,1 (qc,qt,q)>
           *      ...
           */
    PC3   // classical bath, but quantum Qc, Qt, cav-modes & level
          /**
           *
           * Notations:
           *      |g,0,qc,qt (q)>, |e1,0,qc,qt (q)>, |e2,0,qc,qt (q)>
           *      |g,1,qc,qt (q)>, |e1,1,qc,qt (q)>, |e2,1,qc,qt (q)>
           *      ...
           */
};
const std::map<std::string, _enum> _dict = {
    {"pc1", PC1},
    {"pc2", PC2},
    {"pc3", PC3},
};
};  // namespace PyrCavPolicy

class PyrCav_ForceField : public Nad_ForceField {
   public:
    PyrCav_ForceField(const Param& iparm, const int& child);
    PyrCav_ForceField(const Param& iparm);
    PyrCav_ForceField(const std::string& iparm_str) : PyrCav_ForceField(Param::parse(iparm_str)){};

    virtual ~PyrCav_ForceField(){};

    static inline std::string name() { return "pyrincav"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    virtual int ForceField_epes_PC1(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_PC2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_PC3(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

   protected:
    DEFINE_POINTER_PROTECTED(num_real, WCI);
    DEFINE_POINTER_PROTECTED(num_real, ECI);
    DEFINE_POINTER_PROTECTED(num_real, KCI);

    int Fmol, Nmod, Ncav;
    double lcoeff, gcoup;
    double wcav;
    bool first_call = true;

    int fftype;
    std::string ffflag;
};


#endif  // PYRINCAVITY_H

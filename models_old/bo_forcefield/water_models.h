#ifndef Water_ForceField_H
#define Water_ForceField_H

#include <string>

#include "../forcefieldbase.h"

namespace WaterModelPolicy {
enum _enum { qspcfw };
const std::map<std::string, _enum> _dict = {{"#qspcfw", qspcfw}};
}  // namespace WaterModelPolicy

class Water_ForceField : public BO_ForceField {
   public:
    Water_ForceField(const Param& iparm);
    Water_ForceField(const std::string& iparm_str) : Water_ForceField(Param::parse(iparm_str)){};

    virtual ~Water_ForceField(){};

    static inline std::string name() { return "water"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

   protected:
    int Natom;
    int q;
    double boxL, cutoff, ewald_parm;

    DEFINE_POINTER_PROTECTED(num_real, charge_arr);
    DEFINE_POINTER_PROTECTED(num_real, pbox);
    int fftype;
};

#endif
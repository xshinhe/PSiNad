#ifndef SmallMol_ForceField_H
#define SmallMol_ForceField_H

#include "../forcefieldbase.h"

namespace SmallMol_Policy {
enum _enum {
    SM_H2O,
    SM_H2O2,
    SM_NH3,
    SM_paraH2,
    SM_CH2O,
};
const std::map<std::string, int> _dict = {
    {"#h2o", SM_H2O}, {"#h2o2", SM_H2O2}, {"#nh3", SM_NH3}, {"#parah2", SM_paraH2}, {"#ch2o", SM_CH2O},
};

}  // namespace SmallMol_Policy
class SmallMol_ForceField : public BO_ForceField {
   public:
    SmallMol_ForceField(const Param& iparm);
    SmallMol_ForceField(const std::string& iparm_str) : SmallMol_ForceField(Param::parse(iparm_str)){};

    virtual ~SmallMol_ForceField();

    static inline std::string name() { return "smallmol"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

   protected:
    int Natom;
    int Nmole;
    double boxL;
    int fftype;
};

#endif  // SmallMol_ForceField_H

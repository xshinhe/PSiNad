#ifndef MD1D_MODEL_H
#define MD1D_MODEL_H

#include <string>

#include "../forcefieldbase.h"

namespace MD1DPolicy {
enum _enum { MD1D_NONE, MD1D_HO };
const std::map<std::string, _enum> _dict = {
    {"#none", MD1D_NONE},
    {"#ho", MD1D_HO},
};
};  // namespace MD1DPolicy

class MD1D_ForceField : public BO_ForceField {
   public:
    MD1D_ForceField(const Param& iparm, const int& child);
    MD1D_ForceField(const Param& iparm);
    MD1D_ForceField(const std::string& iparm_str) : MD1D_ForceField(Param::parse(iparm_str)){};

    virtual ~MD1D_ForceField();

    static inline std::string name() { return "md1d"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

   protected:
    std::string ffflag;
    int fftype;
};

#endif  // MD1D_MODEL_H

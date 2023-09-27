#ifndef SCTEST_MODEL_H
#define SCTEST_MODEL_H

#include <string>

#include "../forcefieldbase.h"

namespace SCTESTPolicy {
enum _enum { SC1D, SC2D };
const std::map<std::string, _enum> _dict = {
    {"#sc1d", SC1D},
    {"#sc2d", SC2D},
};
};  // namespace SCTESTPolicy

class SCTEST_ForceField : public BO_ForceField {
   public:
    SCTEST_ForceField(const Param& iparm);
    SCTEST_ForceField(const std::string& iparm_str) : SCTEST_ForceField(Param::parse(iparm_str)){};

    virtual ~SCTEST_ForceField();

    static inline std::string name() { return "sctest"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

    int ForceField_npes_SC1D(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                             const int& rdim);

    int ForceField_npes_SC2D(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                             const int& rdim);

   protected:
    std::string ffflag;
    int fftype;
};

#endif  // SCTEST_MODEL_H

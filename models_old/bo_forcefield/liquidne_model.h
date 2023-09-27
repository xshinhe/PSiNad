#ifndef LIQUIDNE_MODEL_H
#define LIQUIDNE_MODEL_H

#include <string>

#include "../forcefieldbase.h"

namespace LiquidneModelPolicy {
enum _enum { liquidne };
const std::map<std::string, _enum> _dict = {{"#liquidne", liquidne}};
}  // namespace LiquidneModelPolicy

class LiquidNe_ForceField : public BO_ForceField {
   public:
    LiquidNe_ForceField(const Param& iparm);

    virtual ~LiquidNe_ForceField();

    static inline std::string name() { return "LiquidNe"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);

   protected:
    int Pnbd, PtN;
};

#endif

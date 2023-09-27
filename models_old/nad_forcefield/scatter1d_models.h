#ifndef Scatter1D_MODELS_H
#define Scatter1D_MODELS_H

#include "../forcefieldbase.h"

namespace Scatter1DPolicy {
enum _enum {
    SAC,   // Tully's Single Avoid Crossing Model
    SAC2,  // Tully's Single Avoid Crossing Model (with slight revision)
    DAC,   // Tully's Double Avoid Crossing Model
    ECR,   // Tully's Extended Coupling Region Model
    DBG,   // (double)
    DAG,   // (double)
    DRN    // (double) ECR
};
const std::map<std::string, _enum> _dict = {
    {"sac", SAC}, {"sac2", SAC2}, {"dac", DAC}, {"ecr", ECR}, {"dbg", DBG}, {"dag", DAG}, {"drn", DRN},
};
};  // namespace Scatter1DPolicy

class Scatter1D_ForceField : public Nad_ForceField {
   public:
    Scatter1D_ForceField(const Param& iparm, const int& child);
    Scatter1D_ForceField(const Param& iparm);
    Scatter1D_ForceField(const std::string& iparm_str) : Scatter1D_ForceField(Param::parse(iparm_str)){};

    virtual ~Scatter1D_ForceField(){};

    static inline std::string name() { return "scatter1d"; }

    virtual int ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                int& eocc, const int& rdim, const int& fdim, const int& icycle);

    virtual int ForceField_spec(num_real* nr, num_real* np, num_real* nm, const int& rdim, const int& fdim);

    virtual int ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P, const int& flag,
                                const int& rdim);


    virtual int Scatter1D_plot(const double& Xrange);

    virtual int ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag, const int& rdim,
                                const int& fdim);

    virtual int ForceField_epes_SAC(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_SAC2(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                     const int& rdim, const int& fdim);

    virtual int ForceField_epes_DAC(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_ECR(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_DBG(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_DAG(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

    virtual int ForceField_epes_DRN(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                    const int& rdim, const int& fdim);

   protected:
    double Xrange;
    int fftype;
    std::string ffflag;
    int spec;
};


#endif  // Scatter1D_MODELS_H

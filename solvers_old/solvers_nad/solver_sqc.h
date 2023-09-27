#ifndef SQC_SOLVER_H
#define SQC_SOLVER_H
#include "nadtraj.h"

#ifndef GAMMA_WIGNER
#define GAMMA_WIGNER(_F) ((std::sqrt((num_real)(_F) + 1) - 1) / (num_real)(_F))
#endif

namespace window {
enum _enum { tri, sqr };
const std::map<std::string, _enum> _dict = {
    {"#tri", tri},
    {"#sqr", sqr},
};
int samp_mvc_window(num_real* mvc, int& iocc, const int& fdim, const int& window_type = window::tri);
int rho_eac(num_complex* rho, num_complex* eac, const int& fdim, const int& window_type = window::tri);
};  // namespace window

class SQC_Solver : public NadTraj_Solver {
   public:
    SQC_Solver(Param iparm, Model* pM);
    SQC_Solver(const std::string& iparm_str, Model* pM) : SQC_Solver(Param::parse(iparm_str), pM){};

    virtual ~SQC_Solver();

    static inline std::string name() { return "sqc"; }

    virtual int init_occ2eac(const int& itraj);
    virtual int init(const int& itraj);
    virtual int kernel0(num_complex* rhox, const int& F);
    virtual int kernelt(num_complex* rhox, const int& F);

   protected:
    int window_type;
    num_real gamma0, gammat;
};

#endif  // SQC_SOLVER_H

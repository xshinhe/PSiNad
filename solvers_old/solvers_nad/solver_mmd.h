#ifndef MMD_SOLVER_H
#define MMD_SOLVER_H
#include "nadtraj.h"

// BETTER_ENUM(phase_enum, int, phase_duni, phase_funi, phase_ddis, phase_fdis)
// BETTER_ENUM(focus_enum, int, focus_mmf, focus_twa, focus_mms1, focus_mms2, focus_mms3)

namespace samp_phase {
enum _enum { duni, funi, ddis, fdis };
const std::map<std::string, _enum> _dict = {
    {"#duni", duni},
    {"#funi", funi},
    {"#ddis", ddis},
    {"#fdis", fdis},
};
};  // namespace samp_phase

namespace focus_action {
enum _enum { mmf, twa, ddd, mms1, mms2, mms3 };
const std::map<std::string, _enum> _dict = {
    {"#mmf", mmf}, {"#twa", twa}, {"#ddd", ddd}, {"#mms1", mms1}, {"#mms2", mms2}, {"#mms3", mms3},
};
};  // namespace focus_action


// enum phase_enum {
//     phase_duni,
//     phase_funi,
//     phase_ddis,
//     phase_fdis,
// };
// const std::map<std::string, phase_enum> phase_dict = {
//     {"#duni", phase_duni},  // dependent uniform
//     {"#funi", phase_funi},  // full uniform
//     {"#ddis", phase_ddis},  // dependent discrete
//     {"#fdis", phase_fdis},  // full discrete
// };
// enum focus_enum {
//     focus_mmf,
//     focus_twa,
//     focus_ddd,
//     focus_mms1,
//     focus_mms2,
//     focus_mms3,
// };
// const std::map<std::string, focus_enum> focus_dict = {
//     {"#mmf", focus_mmf},    // meyer-miller focused
//     {"#twa", focus_twa},    // gdtwa
//     {"#ddd", focus_ddd},    // gdtwa
//     {"#mms1", focus_mms1},  // meyer-miller sampling (uniform)
//     {"#mms2", focus_mms2},  // meyer-miller sampling (linear)
//     {"#mms3", focus_mms3},  // meyer-miller sampling (expotenital)
// };

class MMD_Solver : public NadTraj_Solver {
   public:
    MMD_Solver(Param iparm, Model* pM);
    MMD_Solver(const std::string& iparm_str, Model* pM) : MMD_Solver(Param::parse(iparm_str), pM){};

    virtual ~MMD_Solver();

    static inline std::string name() { return "mmd"; }

    virtual int init(const int& itraj);

    virtual int kernel0(num_complex* rhox, const int& F);

    virtual int kernelt(num_complex* rhox, const int& F);

   protected:
    int phase_type, focus_type;
    int Fref;
    int sums    = 0;
    int occrand = -1;  // random occ

    std::string phase_flag, focus_flag;

    num_real scale, gamma_uu, gamma_ou;
};


#endif  // MMD_SOLVER_H

#ifndef NAD_TCFER_H
#define NAD_TCFER_H

#include "../../utils/definitions.h"

class NAD_TCFer {
   public:
    NAD_TCFer(const int &itype,                // tcf_type
              const int &nspec,                // get specification
              const int &nsamp,                // get nsamp
              const int &N, const int &F = 1,  // problem size
              const int &occ0 = 0              // help to get lentcf
    );

    virtual ~NAD_TCFer(){};

    int Clear();

    int Amount(NAD_TCFer &iT);

    int MPIAmount(NAD_TCFer &iT);

    int Count(const int &isamp, const int &now_spec);

    bool ifrecord(const int &i0, const int &it);

    int report(const std::string &name, const double &unit = 1.0);

    int lentcf;

    bool tcf_reduced;

    DEFINE_POINTER(num_complex, val);
    DEFINE_POINTER(bool, tcf_0t_bool);
    DEFINE_POINTER(num_complex, tcf_0t_val);

   protected:
    std::string header = "";
    int nspec;
    int lsamp;
    int nsamp;
    int ncoll;

    DEFINE_POINTER_PROTECTED(num_complex, coll);
    DEFINE_POINTER_PROTECTED(int, stat);
    DEFINE_POINTER(bool, tcf_0_bool);
    DEFINE_POINTER(bool, tcf_t_bool);
    DEFINE_POINTER(num_complex, tcf_0_val);
    DEFINE_POINTER(num_complex, tcf_t_val);

    int tcf_type;
    std::string tcf_flag;
};

#endif  // NAD_TCFER_H

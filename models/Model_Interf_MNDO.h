#ifndef MODEL_INTERF_MNDO_H
#define MODEL_INTERF_MNDO_H


#include <unistd.h>

#include <cmath>

#include "../core/Kernel.h"

namespace kids {

struct MNDO99KW {
    std::string key;
    int val;
};

using MNDO99KW_map = std::map<std::string, int>;

class Model_Interf_MNDO final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_Interf_MNDO"; }

    Model_Interf_MNDO(){};

    virtual ~Model_Interf_MNDO(){};

   private:
    kids_real temp, beta;
    bool diff_nac;
    std::string init_nuclinp, savename;
    std::vector<std::string> mndo99_data;
    std::string mndo99_keyword, mndo99_comment, mndo99_addition;

    std::vector<MNDO99KW> keyword;  // keyword wrapper

    // integrator
    kids_real *x, *p;

    // model
    int* atoms;
    kids_real* x0;
    kids_real* p0;
    kids_real* x_sigma;
    kids_real* p_sigma;
    kids_real* w;
    kids_real* mass;
    kids_real *vpes, *grad, *hess, *Tmod;
    kids_real *V, *dV;
    kids_real *E, *dE;
    kids_real *nac, *nac_prev;

    int natom;
    int read_flag;
    int ncigrd;
    int iroot;
    bool refer;

    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat);

    int exec_kernel_impl(int stat = -1);

    int parse_mndo99(const std::string& mndo99inp);
    std::string new_keyword(const MNDO99KW_map& newkeyword);
    int new_task(const std::string& file, const std::string& task_flag);
    int track_nac_sign();
    int parse_standard(const std::string& log);
    int parse_hessian(const std::string& log);
    int calc_normalmode();
    int calc_samp();
    int calc_scan();
};
};  // namespace kids

#endif  // MODEL_INTERF_MNDO_H
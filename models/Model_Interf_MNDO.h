#ifndef MODEL_INTERF_MNDO_H
#define MODEL_INTERF_MNDO_H


#include <unistd.h>

#include <cmath>

#include "../core/Kernel.h"

namespace PROJECT_NS {

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
    num_real temp, beta;
    bool diff_nac;
    std::string init_nuclinp, savename;
    std::vector<std::string> mndo99_data;
    std::string mndo99_keyword, mndo99_comment, mndo99_addition;

    std::vector<MNDO99KW> keyword;  // keyword wrapper

    // integrator
    num_real *x, *p;

    // model
    int* atoms;
    num_real* x0;
    num_real* p0;
    num_real* x_sigma;
    num_real* p_sigma;
    num_real* w;
    num_real* mass;
    num_real *vpes, *grad, *hess, *Tmod;
    num_real *V, *dV;
    num_real *E, *dE;
    num_real *nac, *nac_prev;

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
};  // namespace PROJECT_NS

#endif  // MODEL_INTERF_MNDO_H
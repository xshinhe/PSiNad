#ifndef MODEL_INTERF_MNDO_H
#define MODEL_INTERF_MNDO_H


#include <unistd.h>

#include <cmath>

#include "../core/Kernel.h"

namespace PROJECT_NS {

struct MNDOKW {
    std::string key;
    int val;
};

using MNDOKW_map = std::map<std::string, int>;

class Model_Interf_MNDO final : public Kernel {
   public:
    inline virtual const std::string name() { return "Model_Interf_MNDO"; }

    Model_Interf_MNDO(){};

    virtual ~Model_Interf_MNDO(){};

   private:
    kids_real temp, beta;
    bool diff_nac;

    std::string exec_file;
    std::string init_nuclinp;
    std::string savename;

    std::vector<std::string> mndo_data;
    std::string mndo_keyword;
    std::string mndo_comment;
    std::string mndo_addition;

    std::vector<MNDOKW> keyword;  // keyword wrapper

    std::string task_control;
    std::string directory;

    bool classical_bath;

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
    kids_real *T, *E, *dE;
    kids_real *nac, *nac_prev;

    kids_real *f_r, *f_p, *f_rp;

    int natom;
    int read_flag;
    int nciref;
    int ncigrd;
    int iroot;
    int lroot;
    bool refer;
    bool* frez_ptr;

    void read_param_impl(Param* PM);

    void init_data_impl(DataSet* DS);

    void init_calc_impl(int stat);

    int exec_kernel_impl(int stat = -1);

    int parse_mndo(const std::string& mndoinp);
    std::string new_keyword(const MNDOKW_map& newkeyword);
    int new_task(const std::string& file, const std::string& task_flag);
    int track_nac_sign();
    int parse_standard(const std::string& log, int stat);
    int parse_hessian(const std::string& log);
    int parse_hessian2(const std::string& log);
    int calc_normalmode();
    int calc_samp();
    int calc_scan();
};
};  // namespace PROJECT_NS

#endif  // MODEL_INTERF_MNDO_H
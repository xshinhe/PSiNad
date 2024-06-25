#ifndef MODEL_INTERF_MNDO_H
#define MODEL_INTERF_MNDO_H


#include <unistd.h>

#include <cmath>

#include "kids/Model.h"

namespace PROJECT_NS {

struct MNDOKW {
    std::string key;
    std::string val;
};

using MNDOKW_map = std::map<std::string, std::string>;

class Model_Interf_MNDO final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_Interf_MNDO(){};

    virtual ~Model_Interf_MNDO(){};

   private:
    kids_real temp, beta;
    bool      diff_nac;

    std::string exec_file;
    std::string init_nuclinp;
    std::string savename;

    std::vector<std::string> mndo_data;
    std::string              mndo_keyword;
    std::string              mndo_comment;
    std::string              mndo_addition;

    std::vector<MNDOKW> keyword;  // keyword wrapper

    std::string task_control;
    std::string directory;

    bool classical_bath;

    // integrator
    kids_real *x, *p;

    // model
    int*       atoms;
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

    int  natom;
    int  read_flag;
    int  nciref;
    int  ncigrd;
    int  iroot;
    int  lroot;
    bool refer;

    kids_bint* succ_ptr;
    kids_bint* frez_ptr;
    kids_bint* last_attempt_ptr;
    int*       fail_type_ptr;

    void setInputParam_impl(std::shared_ptr<Param> PM);

    void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    Status& initializeKernel_impl(Status& stat);

    Status& executeKernel_impl(Status& stat);

    int         parse_mndo(const std::string& mndoinp);
    std::string new_keyword(const MNDOKW_map& newkeyword);
    int         new_task(const std::string& file, const std::string& task_flag);
    int         track_nac_sign();
    Status&     parse_standard(const std::string& log, Status& stat);
    int         parse_hessian(const std::string& log);
    int         parse_hessian2(const std::string& log);
    int         calc_normalmode();
    int         calc_samp();
    int         calc_scan();
};
};  // namespace PROJECT_NS

#endif  // MODEL_INTERF_MNDO_H
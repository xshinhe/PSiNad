#ifndef MODEL_QMInterface_H
#define MODEL_QMInterface_H

#include "kids/Model.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(QMPolicy,  //
              GAUSSIAN,  //
              ORCA,      //
              MNDO,      //
              BAGEL,     //
              MOLCAS,    //
              NONE);     //

class Model_QMInterface final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_QMInterface(){};

    virtual ~Model_QMInterface(){};

   private:
    QMPolicy::_type qm_type;

    std::string pykids_path;
    std::string qm_config_in;
    std::string config_content;
    std::string exec_file;

    kids_real   temp, beta;
    kids_real   ener_refered;
    kids_real   time_unit;
    bool        diff_nac;
    std::string init_nuclinp;
    std::string savename;
    std::string task_control;

    bool save_every_calc;
    bool save_every_step;

    // integrator
    kids_real *x, *p;
    int*       atoms;
    kids_real* x0;
    kids_real* p0;
    kids_real* x_sigma;
    kids_real* p_sigma;
    kids_real* w;
    kids_real* mass;
    kids_real *vpes, *grad, *hess, *Tmod;
    kids_real *V, *dV;
    kids_real *T, *eig, *dE;
    kids_real *nac, *nac_prev;
    kids_real *f_r, *f_p, *f_rp;
    kids_real *dt_ptr, *t_ptr;
    kids_int*  istep_ptr;

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

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);

    int track_nac_sign();
};

};      // namespace PROJECT_NS
#endif  // MODEL_QMInterface_H
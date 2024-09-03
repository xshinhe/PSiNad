#ifndef MODEL_QMMMInterface_H
#define MODEL_QMMMInterface_H

#include "kids/Model.h"
#include "kids/Policy.h"

namespace PROJECT_NS {

DEFINE_POLICY(QMMMPolicy,  //
              GAUSSIAN,    //
              ORCA,        //
              MNDO,        //
              BAGEL,       //
              MOLCAS,      //
              NONE);       //

class Model_QMMMInterface final : public Model {
   public:
    virtual const std::string getName();

    virtual int getType() const;

    Model_QMMMInterface(){};

    virtual ~Model_QMMMInterface(){};

   private:
    QMMMPolicy::_type QMMM_type;

    std::string kidsqmmm_path;
    std::string qmmm_config_in;
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
    int  sstep_dataset;

    // integrator
    span<kids_real> x, p;
    span<kids_int>  atoms;
    span<kids_real> x0;
    span<kids_real> p0;
    span<kids_real> x_sigma;
    span<kids_real> p_sigma;
    span<kids_real> w;
    span<kids_real> mass;
    span<kids_real> vpes, grad, hess, Tmod;
    span<kids_real> V, dV;
    span<kids_real> T, eig, dE;
    span<kids_real> nac, nac_prev;
    span<kids_real> f_r, f_p, f_rp;
    span<kids_real> dt_ptr, t_ptr;
    span<kids_int>  istep_ptr;

    int  natom;
    int  read_flag;
    int  nciref;
    int  ncigrd;
    int  iroot;
    int  lroot;
    bool refer;

    int try_level;

    span<kids_bint> succ_ptr;
    span<kids_bint> frez_ptr;
    span<kids_bint> last_attempt_ptr;
    span<kids_int>  fail_type_ptr;

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);

    int track_nac_sign();
};

};      // namespace PROJECT_NS
#endif  // MODEL_QMMMInterface_H
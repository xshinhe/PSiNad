#ifndef MODEL_QMInterface_H
#define MODEL_QMInterface_H

#include "psnd/Model.h"
#include "psnd/Policy.h"

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

    Model_QMInterface() {};

    virtual ~Model_QMInterface() {};

   private:
    QMPolicy::_type qm_type;

    std::string pypsnd_path;
    std::string qm_string;
    std::string qm_config_in;
    std::string config_content;
    std::string exec_file;

    psnd_real   temp, beta;
    psnd_real   ener_refered;
    psnd_real   time_unit;
    bool        diff_nac;
    std::string init_nuclinp;
    std::string savename;
    std::string task_control;

    bool save_every_calc;
    bool save_every_step;
    int  sstep_dataset;
    bool use_state_detection;

    // integrator
    span<psnd_real> x, p;
    span<psnd_int>  atoms;
    span<psnd_real> x0;
    span<psnd_real> p0;
    span<psnd_real> x_sigma;
    span<psnd_real> p_sigma;
    span<psnd_real> w;
    span<psnd_real> mass;
    span<psnd_real> vpes, grad, hess, Tmod;
    span<psnd_real> V, dV;
    span<psnd_real> T, eig, dE;
    span<psnd_real> nac, nac_prev;
    span<psnd_real> f_r, f_p, f_rp;
    span<psnd_real> dt_ptr, t_ptr;
    span<psnd_int>  istep_ptr;

    int  natom;
    int  read_flag;
    int  nciref;
    int  ncigrd;
    int  iroot;
    int  lroot;
    bool refer;

    int try_level;

    span<psnd_bint> succ_ptr;
    span<psnd_bint> frez_ptr;
    span<psnd_bint> last_attempt_ptr;
    span<psnd_int>  fail_type_ptr;

    virtual void    setInputParam_impl(std::shared_ptr<Param> PM);
    virtual void    setInputDataSet_impl(std::shared_ptr<DataSet> DS);
    virtual Status& initializeKernel_impl(Status& stat);
    virtual Status& executeKernel_impl(Status& stat);

    int track_nac_sign();
};

};  // namespace PROJECT_NS
#endif  // MODEL_QMInterface_H

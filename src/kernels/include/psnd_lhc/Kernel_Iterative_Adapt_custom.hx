/**@file        Kernel_Iterative_Adapt_custom.h
 * @brief       custom copy of Kernel_Iterative_Adapt that supports
 *              reducing timestep when the minimum electronic energy gap
 *              is below a user-set threshold.
 *
 * Note: this is a user-level copy so original sources are unchanged.
 */

#ifndef Kernel_Iterative_Adapt_custom_H
#define Kernel_Iterative_Adapt_custom_H

#include <chrono>

#include "psnd/Kernel.h"

namespace PROJECT_NS {

class Kernel_Iterative_Adapt_custom final : public Kernel {
   public:
    Kernel_Iterative_Adapt_custom() { enable_call_child = false; }

    virtual const std::string getName();

    virtual int getType() const;

   private:
    int             nstep, sstep, nsamp;
    double          t0, tend, dt0;
    span<psnd_real> t, dt;
    span<psnd_bint> at_condition;

    int            msize;
    span<psnd_int> tsize, dtsize, last_tried_dtsize;
    span<psnd_int> istep, isamp;
    int            nbackup;
    double         time_unit;

    bool                                               restart;
    bool                                               use_exchange;
    int                                                exchange_root;
    int                                                exchange_num;
    double                                             exchange_time;
    std::chrono::time_point<std::chrono::steady_clock> ex_begin;

    // Adaptive dt parameters (read from Param)
    double adapt_gap_threshold = 0.0;  // if <= 0, feature disabled
    double adapt_dt_factor     = 0.5;  // multiply dt by this factor when gap small

    // pointer to electronic eigenvalues (if present)
    span<psnd_real> eig;   // DATA::model::rep::eig
    span<psnd_real> dE;    // DATA::model::rep::dE (per-state energy differences)

    const std::vector<std::string> backup_fields = {"x", "p", "U", "occ_nuc", "f", "Ekin", "Epot", "rho_ele", "c"};

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual Status& initializeKernel_impl(Status& stat);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status& executeKernel_impl(Status& stat);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Iterative_Adapt_custom_H

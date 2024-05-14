#ifndef Handler_H
#define Handler_H

#include "kids/Kernel.h"
#include "kids/Param.h"

namespace PROJECT_NS {

/**
 * this class provides a control of simulation
 */
class Handler final {
   public:
    Handler(const std::string& solver_name, const std::string& model_name);

    virtual ~Handler(){};

    int run(std::shared_ptr<Param>& PM);

    int run_single(std::shared_ptr<Param>& PM);

    int run_single_mpi(std::shared_ptr<Param>& PM);

    int run_parallel(std::shared_ptr<Param>& PM);

    int run_sampling(std::shared_ptr<Param>& PM);  //{ return 0; }

    int run_help(std::shared_ptr<Param>& PM) { return 0; }

    int run_help_param(std::shared_ptr<Param>& PM) { return 0; }

    int run_help_dataset(std::shared_ptr<Param>& PM);  // { return 0; }

   private:
    std::shared_ptr<Kernel> model;
    std::shared_ptr<Kernel> solver;
};

};  // namespace PROJECT_NS

#endif  // Handler_H
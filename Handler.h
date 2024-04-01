#ifndef Handler_H
#define Handler_H

#include "core/Kernel.h"
#include "core/Param.h"

namespace PROJECT_NS {

/**
 * this class provides a control of simulation
 */
class Handler final {
   public:
    Handler(const std::string& solver_name, const std::string& model_name);

    virtual ~Handler(){};

    int run(Param* PM);

    int run_single(Param* PM);

    int run_single_mpi(Param* PM);

    int run_parallel(Param* PM);

    int run_sampling(Param* PM);  //{ return 0; }

    int run_help(Param* PM) { return 0; }

    int run_help_param(Param* PM) { return 0; }

    int run_help_dataset(Param* PM);  // { return 0; }

   private:
    std::shared_ptr<Kernel> model;
    std::shared_ptr<Kernel> solver;
};

};  // namespace PROJECT_NS

#endif  // Handler_H
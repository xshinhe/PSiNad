#ifndef Handler_H
#define Handler_H

#include "core/Kernel.h"
#include "core/Param.h"

namespace PROJECT_NS {


class Handler final {
   public:
    // enum hdlr_t {
    //     parallel,
    //     single,
    //     single_mpi,
    //     sampling,
    //     help,
    //     help_param,
    //     help_dataset,
    // };

    // static const std::map<std::string, hdlr_t> _dict = {
    //     {"parallel", parallel},
    //     {"single", single},
    //     {"single_mpi", single_mpi},
    //     {"sampling", sampling},
    // };

    Handler(const std::string& solver_name, const std::string& model_name);

    virtual ~Handler(){};

    int run(Param* PM);

    int run_single(Param* PM);

    int run_single_mpi(Param* PM);

    int run_parallel(Param* PM);

    int run_sampling(Param* PM);  //{ return 0; }

    int run_help(Param* PM) { return 0; }

    int run_help_param(Param* PM) { return 0; }

    int run_help_dataset(Param* PM) { return 0; }

   private:
    std::shared_ptr<Kernel> model;
    std::shared_ptr<Kernel> solver;
};

};  // namespace PROJECT_NS

#endif  // Handler_H
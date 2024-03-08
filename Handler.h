#ifndef Handler_H
#define Handler_H

#include "core/Kernel.h"
#include "core/Param.h"

namespace PROJECT_NS {


class Handler final {
   public:
    enum hdlr_t {
        single,
        multiple,
    };

    Handler(hdlr_t itype, const std::string& solver_name, const std::string& model_name);

    virtual ~Handler(){};

    int run(Param* P);

    int run_single(Param* P);

    int run_multiple(Param* P);

   private:
    hdlr_t type = hdlr_t::single;
    std::shared_ptr<Kernel> model;
    std::shared_ptr<Kernel> solver;
};

};  // namespace PROJECT_NS

#endif  // Handler_H
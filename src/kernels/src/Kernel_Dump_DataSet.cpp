#include "kids/Kernel_Dump_DataSet.h"

#include "kids/debug_utils.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

const std::string Kernel_Dump_DataSet::getName() { return "Kernel_Dump_DataSet"; }

int Kernel_Dump_DataSet::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Dump_DataSet::setInputParam_impl(std::shared_ptr<Param> PM) {
    fn             = _param->get_string({"solver.dump", "dump"}, LOC(), "final");
    hdlr_str       = _param->get_string({"solver.handler", "handler"}, LOC(), "");
    dump_only_init = _param->get_bool({"solver.dump_only_init"}, LOC(), false);
}

Status& Kernel_Dump_DataSet::initializeKernel_impl(Status& stat) { return stat; }

Status& Kernel_Dump_DataSet::executeKernel_impl(Status& stat) {
    if (fn == "" || fn == "null") return stat;
    try {
        std::ofstream ofs{utils::concat(directory, "/", fn, stat.icalc, ".ds")};
        if (dump_only_init) {
            _dataset->dump_match(ofs, "init");
        } else {
            _dataset->dump(ofs);
        }
        ofs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return stat;
}

};  // namespace PROJECT_NS

#include "kids/Kernel_Dump_DataSet.h"

#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

const std::string Kernel_Dump_DataSet::getName() { return "Kernel_Dump_DataSet"; }

int Kernel_Dump_DataSet::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Dump_DataSet::setInputParam_impl(std::shared_ptr<Param> PM) {
    fn       = _param->get_string({"solver.dump", "dump"}, LOC(), "final");
    hdlr_str = _param->get_string({"solver.handler", "handler"}, LOC(), "");
}

Status& Kernel_Dump_DataSet::initializeKernel_impl(Status& stat) {
    if (hdlr_str == "sampling") {
        try {
            std::ofstream ofs{utils::concat(directory, "/samp", stat.icalc, ".ds")};
            _dataset->dump(ofs);
            ofs.close();
        } catch (std::runtime_error& e) { throw kids_error(fn); }
    }
    return stat;
}

Status& Kernel_Dump_DataSet::executeKernel_impl(Status& stat) {
    if (fn == "" || fn == "null") return stat;
    try {
        std::ofstream ofs{utils::concat(directory, "/", fn, stat.icalc, ".ds")};
        _dataset->dump(ofs);
        ofs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return stat;
}

};  // namespace PROJECT_NS

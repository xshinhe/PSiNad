#include "kids/Kernel_Dump_DataSet.h"

#include "kids/linalg.h"

namespace PROJECT_NS {

void Kernel_Dump_DataSet::setInputParam_impl(std::shared_ptr<Param>& PM) {
    fn        = PM->get_string("dump", LOC(), "final");
    hdlr_str  = PM->get_string("handler", LOC());
    directory = PM->get_string("directory", LOC());
    assert(directory != "");
}

Status& Kernel_Dump_DataSet::initializeKernel_impl(Status& stat) {
    if (hdlr_str == "sampling") {
        try {
            std::ofstream ofs{utils::concat(directory, "/samp", stat, ".ds")};
            _dataset->dump(ofs);
            ofs.close();
        } catch (std::runtime_error& e) { throw kids_error(fn); }
    }
}

Status& Kernel_Dump_DataSet::executeKernel_impl(Status& stat) {
    if (fn == "" || fn == "null") return 0;
    try {
        std::ofstream ofs{utils::concat(directory, "/", fn, stat, ".ds")};
        _dataset->dump(ofs);
        ofs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return 0;
}

};  // namespace PROJECT_NS

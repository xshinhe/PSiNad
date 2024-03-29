#include "Kernel_Dump_DataSet.h"

#include "../core/linalg.h"

namespace PROJECT_NS {

void Kernel_Dump_DataSet::read_param_impl(Param* PM) {
    fn        = PM->get<std::string>("dump", LOC(), "final");
    hdlr_str  = PM->get<std::string>("handler", LOC());
    directory = PM->get<std::string>("directory", LOC());
    assert(directory != "");
}

void Kernel_Dump_DataSet::init_calc_impl(int stat) {
    if (hdlr_str == "sampling") {
        try {
            std::ofstream ofs{utils::concat(directory, "/samp", stat, ".ds")};
            _DataSet->dump(ofs);
            ofs.close();
        } catch (std::runtime_error& e) { throw kids_error(fn); }
    }
}

int Kernel_Dump_DataSet::exec_kernel_impl(int stat) {
    if (fn == "" || fn == "null") return 0;
    try {
        std::ofstream ofs{utils::concat(directory, "/", fn, stat, ".ds")};
        _DataSet->dump(ofs);
        ofs.close();
    } catch (std::runtime_error& e) { throw kids_error(fn); }
    return 0;
}

};  // namespace PROJECT_NS

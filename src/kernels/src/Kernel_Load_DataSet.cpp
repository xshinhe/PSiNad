#include "psnd/Kernel_Load_DataSet.h"

#include "psnd/hash_fnv1a.h"
#include "psnd/linalg.h"
#include "psnd/macro_utils.h"
#include "psnd/vars_list.h"

namespace PROJECT_NS {

const std::string Kernel_Load_DataSet::getName() { return "Kernel_Load_DataSet"; }

int Kernel_Load_DataSet::getType() const { return utils::hash(FUNCTION_NAME); }

void Kernel_Load_DataSet::setInputParam_impl(std::shared_ptr<Param> PM) {
    load_fn = _param->get_string({"load", "solver.load"}, LOC(), "");
}

Status& Kernel_Load_DataSet::initializeKernel_impl(Status& stat) {
    if (load_fn == "") return stat;
    try {
        auto        ipos         = load_fn.find(":");
        std::string load_fn_file = (ipos == std::string::npos) ? load_fn : load_fn.substr(0, ipos);
        if (load_fn_file.find(".ds") != std::string::npos) {
            std::ifstream ifs{load_fn_file};
            _dataset_load = std::shared_ptr<DataSet>(new DataSet());
            _dataset_load->load(ifs);
            ifs.close();
        } else {
            std::ifstream ifs{utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")};
            _dataset_load->load(ifs);
            ifs.close();
        }
        syncDataSetLoad(_dataset_load);
        if (load_fn.find(":continue") != std::string::npos) {
            if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));
            std::istringstream iss(_dataset_load->repr());
            _dataset->load(iss);
        } else if (load_fn.find(":restart") != std::string::npos) {
            if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));
            std::istringstream iss(_dataset_load->repr());
            int                nsamp = _dataset->def(VARIABLE<psnd_int>("control.nsamp", &Dimension::shape_1, "@"))[0];
            _dataset->load_reframe(iss, nsamp);
        }
    } catch (std::runtime_error& e) { throw psnd_error(load_fn); }
    return stat;
}

};  // namespace PROJECT_NS

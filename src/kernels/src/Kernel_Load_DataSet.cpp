#include "kids/Kernel_Load_DataSet.h"

#include "kids/hash_fnv1a.h"
#include "kids/linalg.h"
#include "kids/macro_utils.h"
#include "kids/vars_list.h"
#include <filesystem>

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
            
            std::cout << "Loading DataSet from file: " << load_fn_file << "\n";
            // 先确认有没有这个文件
            if (!std::filesystem::exists(load_fn_file)) {
                throw kids_error(utils::concat("File does not exist: ", load_fn_file));
            }
            std::ifstream ifs{load_fn_file};
            
            if (!ifs.is_open()) {
                throw kids_error(utils::concat("Cannot open file: ", load_fn_file));
            }
            _dataset_load = std::shared_ptr<DataSet>(new DataSet());
            _dataset_load->load(ifs);
            ifs.close();
        } else {

            std::cout << "Loading DataSet from directory: " << utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds") << "\n";

            // 先确认有没有这个文件
            if (!std::filesystem::exists(utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds"))) {
                throw kids_error(utils::concat("File does not exist: ", utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")));
            }

            std::ifstream ifs{utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")};
            if (!ifs.is_open()) {
                throw kids_error(utils::concat("Cannot open file: ", utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")));
            }
            _dataset_load->load(ifs);
            ifs.close();
        }
        syncDataSetLoad(_dataset_load);
        if (load_fn.find(":continue") != std::string::npos) {
            if (_dataset_load == nullptr) throw kids_error(utils::concat(LOC(), ": DataSet Load error"));
            std::istringstream iss(_dataset_load->repr());
            _dataset->load(iss);
        } else if (load_fn.find(":restart") != std::string::npos) {
            if (_dataset_load == nullptr) throw kids_error(utils::concat(LOC(), ": DataSet Load error"));
            std::istringstream iss(_dataset_load->repr());
            int                nsamp = _dataset->def(VARIABLE<kids_int>("control.nsamp", &Dimension::shape_1, "@"))[0];
            _dataset->load_reframe(iss, nsamp);
        }
    } catch (std::runtime_error& e) { throw kids_error(load_fn); }
    return stat;
}

};  // namespace PROJECT_NS

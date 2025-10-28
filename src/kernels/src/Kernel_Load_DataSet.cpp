#include "psnd/Kernel_Load_DataSet.h"

// @blame thirdpart/filesystem is a wrapper of <filesystem> and compile with previous version of c++ (c++11)
#include <filesystem>

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
            std::cout << "Loading DataSet from file: " << load_fn_file << "\n";
            // 先确认有没有这个文件 @ you can use ghc::filesystem instead for std=c++11
            if (!std::filesystem::exists(load_fn_file)) {
                throw psnd_error(utils::concat("File does not exist: ", load_fn_file));
            }

            if (load_fn_file.find(".bin.ds") != std::string::npos) {
                std::cout << "Loading dataset in binary format" << "\n";
                _dataset_load = std::shared_ptr<DataSet>(new DataSet());
                _dataset_load->load_binary(load_fn_file);
            } else {
                std::ifstream ifs{load_fn_file};

                if (!ifs.is_open()) { throw psnd_error(utils::concat("Cannot open file: ", load_fn_file)); }
                _dataset_load = std::shared_ptr<DataSet>(new DataSet());
                _dataset_load->load(ifs);
                // std::cout << _dataset_load->repr() << "\n";
                ifs.close();
            }

        } else {
            std::cout << "Loading DataSet from directory: "
                      << utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds") << "\n";

            // 先确认有没有这个文件
            if (!std::filesystem::exists(utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds"))) {
                throw std::runtime_error(utils::concat("File does not exist: ",
                                                       utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")));
            }

            std::ifstream ifs{utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")};
            if (!ifs.is_open()) {
                throw psnd_error(utils::concat("Cannot open file: ",
                                               utils::concat(directory, "/", load_fn_file, stat.icalc, ".ds")));
            }

            _dataset_load = std::shared_ptr<DataSet>(new DataSet());
            _dataset_load->load(ifs);

            ifs.close();
        }
        std::cout << "DataSet loaded successfully.\n";
        syncDataSetLoad(_dataset_load);
        std::cout << "DataSet synchronized successfully.\n";
        if (load_fn.find(":continue") != std::string::npos) {
            if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));

            std::cout << "Continuing DataSet from loaded state...\n";
            std::istringstream iss(_dataset_load->repr());
            _dataset->load(iss);
        } else if (load_fn.find(":restart") != std::string::npos) {
            if (_dataset_load == nullptr) throw psnd_error(utils::concat(LOC(), ": DataSet Load error"));

            std::cout << "Reframing DataSet based on nsamp...\n";
            std::istringstream iss(_dataset_load->repr());
            int                nsamp = _dataset->def(VARIABLE<psnd_int>("control.nsamp", &Dimension::shape_1, "@"))[0];
            std::cout << "nsamp = " << nsamp << "\n";
            _dataset->load_reframe(iss, nsamp);
        }
    } catch (const psnd_error& e) {
        // 重新抛出 psnd_error，保持友好的错误信息
        throw;
    } catch (const std::filesystem::filesystem_error& e) {
        throw psnd_error(
            utils::concat("❌ 文件系统错误: ", e.what(),
                          "\n请检查:\n  1. 文件路径是否正确\n  2. 文件权限是否足够\n  3. 磁盘空间是否充足"));
    } catch (const std::ios_base::failure& e) {
        throw psnd_error(
            utils::concat("❌ 文件读写错误: ", e.what(),
                          "\n请检查:\n  1. 文件是否存在且可读\n  2. 文件是否被其他程序占用\n  3. 磁盘是否有足够空间"));
    } catch (const std::runtime_error& e) {
        std::string original_error = e.what();
        throw psnd_error(
            utils::concat("❌ 数据集加载失败: ", original_error,
                          "\n这可能是由于:\n  1. 重启文件格式错误\n  2. 数据不完整或损坏\n  3. 配置参数不匹配"));
    } catch (const std::exception& e) {
        throw psnd_error(utils::concat("❌ 未知错误: ", e.what(),
            "\n如果问题持续出现，请联系开发团队并提供完整的错误信息"));
    } catch (...) {
        throw psnd_error("❌ 未知严重错误发生，程序终止。\n请联系开发团队并提供尽可能多的上下文信息以协助排查问题。");
    }
    return stat;
}

};  // namespace PROJECT_NS

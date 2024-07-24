#include "kids/Kernel.h"

#include <chrono>
#include <ctime>

#include "kids/Exception.h"
#include "kids/hash_fnv1a.h"
#include "kids/macro_utils.h"

namespace PROJECT_NS {

Kernel::Kernel(const std::string& customized_name) : kernel_name{customized_name} {
    kernel_id      = Kernel::getKernels().size();
    kernel_type    = getType();
    max_align_size = getName().size();
    Kernel::getKernels().push_back(this);
};

Kernel::~Kernel(){};

void Kernel::setTiming(bool is_timing_in) {
    is_timing = is_timing_in;
    for (auto& pkernel : _child_kernels) pkernel->setTiming(is_timing_in);
}

void Kernel::setRuleSet(std::shared_ptr<RuleSet> ruleset) {
    _ruleset = ruleset;
    for (auto& pkernel : _child_kernels) pkernel->setRuleSet(ruleset);
}

void Kernel::setInputParam(std::shared_ptr<Param> PM) {
    _param     = PM;
    is_timing  = _param->get_bool({"solver.timing", "timing"}, LOC(), is_timing);
    montecarlo = _param->get_int({"solver.M"}, LOC(), 1);
    directory  = _param->get_string({"solver.directory", "directory"}, LOC(), "default");
    setInputParam_impl(PM);
    for (auto& pkernel : _child_kernels) pkernel->setInputParam(PM);
}

void Kernel::setInputDataSet(std::shared_ptr<DataSet> DS) {
    if (!_param) throw kids_error("Param must be passed before");
    _dataset = DS;
    setInputDataSet_impl(DS);
    for (auto& pkernel : _child_kernels) pkernel->setInputDataSet(DS);
}

std::shared_ptr<Param> Kernel::getParam() const { return _param; }

std::shared_ptr<DataSet> Kernel::getDataSet() const { return _dataset; }

Status& Kernel::initializeKernel(Status& stat) {
    if (!_dataset) throw kids_error("DataSet must be passed before");
    // std::cout << "init: " << LOC() << getName() << "\n";
    // @todo: Consider if the load option is available and ensure it is not overwritten by this function.
    initializeKernel_impl(stat);
    for (auto& pkernel : _child_kernels) pkernel->initializeKernel(stat);
    count_calc++;  // random only called once!!
    return stat;
}

Status& Kernel::executeKernel(Status& stat) {
    if (!_dataset) throw kids_error("DataSet must be passed before");
    // if (!_ruleset & !has_parent) std::cerr << "run without rules\n";

    // std::cout << "exec: " << LOC() << getName() << "\n";
    std::chrono::time_point<std::chrono::steady_clock> begin, end;
    if (is_timing) begin = std::chrono::steady_clock::now();
    {
        executeKernel_impl(stat);
        if (enable_call_child)
            for (auto& pkernel : _child_kernels) pkernel->executeKernel(stat);
    }
    if (is_timing) {
        end = std::chrono::steady_clock::now();
        exec_time += static_cast<std::chrono::duration<double>>(end - begin).count();
    }
    count_exec++;
    return stat;
}

Status& Kernel::finalizeKernel(Status& stat) {
    finalizeKernel_impl(stat);
    for (auto& pkernel : _child_kernels) pkernel->finalizeKernel(stat);
    return stat;
}

int Kernel::getType() const { return utils::hash(FUNCTION_NAME); }

int Kernel::getID() const { return kernel_id; }

bool Kernel::operator==(const Kernel& ker) { return getID() == ker.getID(); }

Kernel& Kernel::appendChild(std::shared_ptr<Kernel> ker) {
    ker->_order_in_parent = _child_kernels.size();
    _child_kernels.push_back(ker);
    ker->_parent_kernel = this;
    ker->has_parent     = true;

    if (ker->_ruleset) {
        _ruleset = ker->_ruleset;
        Kernel* current;
        size_t  _dummy;
        std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
        while (current != nullptr) {
            current->_ruleset         = ker->_ruleset;
            std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
        }
    }

    // shouldn't be here!!
    depth = std::max(depth, ker->depth + 1);
    if (ker->getName().size() > max_align_size) { max_align_size = ker->getName().size(); }

    return *this;
}

Kernel& Kernel::insertAt(std::vector<std::size_t> indexes, std::shared_ptr<Kernel> ker) {
    if (indexes.size() <= 0) throw kids_error("Bad indexes for accessing Kernel");
    std::size_t idx0 = indexes[0];
    if (idx0 > _child_kernels.size()) throw kids_error("Bad indexes for accessing Kernel");
    if (indexes.size() == 1) {
        ker->_order_in_parent = idx0;
        _child_kernels.insert(_child_kernels.begin() + idx0, ker);
        ker->_parent_kernel = this;
        ker->has_parent     = true;
        if (ker->_ruleset) {
            _ruleset = ker->_ruleset;
            Kernel* current;
            size_t  _dummy;
            std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
            while (current != nullptr) {
                current->_ruleset         = ker->_ruleset;
                std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
            }
        }
    } else {
        indexes.erase(indexes.begin());
        _child_kernels[idx0]->insertAt(indexes, ker);
    }
    return *this;
}

Kernel& Kernel::removeAt(std::vector<std::size_t> indexes) {
    if (indexes.size() <= 0) throw kids_error("Bad indexes for accessing Kernel");
    std::size_t idx0 = indexes[0];
    if (idx0 >= _child_kernels.size()) throw kids_error("Bad indexes for accessing Kernel");
    if (indexes.size() == 1) {
        auto& ker           = _child_kernels[idx0];
        ker->_parent_kernel = nullptr;
        ker->has_parent     = false;
        _child_kernels.erase(_child_kernels.begin() + idx0);
    } else {
        indexes.erase(indexes.begin());
        _child_kernels[idx0]->removeAt(indexes);
    }
    return *this;
}

/**
 * @brief build tree structure of the kernel
 */
Kernel& Kernel::updateAt(std::vector<std::size_t> indexes, std::shared_ptr<Kernel> ker) {
    if (indexes.size() <= 0) throw kids_error("Bad indexes for accessing Kernel");
    std::size_t idx0 = indexes[0];
    if (idx0 >= _child_kernels.size()) throw kids_error("Bad indexes for accessing Kernel");
    if (indexes.size() == 1) {
        auto& kerold           = _child_kernels[idx0];
        kerold->_parent_kernel = nullptr;
        kerold->has_parent     = false;

        ker->_order_in_parent = idx0;
        _child_kernels[idx0]  = ker;
        ker->_parent_kernel   = this;
        ker->has_parent       = true;

        if (ker->_ruleset) {
            _ruleset = ker->_ruleset;
            Kernel* current;
            size_t  _dummy;
            std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
            while (current != nullptr) {
                current->_ruleset         = ker->_ruleset;
                std::tie(current, _dummy) = getLastParentKernelAndChildOrder();
            }
        }
    } else {
        indexes.erase(indexes.begin());
        _child_kernels[idx0]->updateAt(indexes, ker);
    }
    return *this;
}

std::tuple<Kernel*, std::size_t> Kernel::getLastParentKernelAndChildOrder() {
    if (has_parent) return std::make_tuple(_parent_kernel, _order_in_parent);
    return std::make_tuple(nullptr, 0);
}

std::shared_ptr<RuleSet> Kernel::getRuleSet() { return _ruleset; }

const std::string Kernel::serializeKernel(const Kernel& ker) { return ""; }

std::shared_ptr<Kernel> Kernel::deserializeKernel(const std::string& str) { return nullptr; }

const std::string Kernel::generateInformationString(double total_time, int current_layer, int total_depth,
                                                    int total_align_size) {
    std::stringstream ss;

    if (total_depth == 0) total_depth = depth;
    if (total_align_size == 0) total_align_size = max_align_size;
    if (total_time < 0) total_time = exec_time;

    if (current_layer == 0) {
        ss << std::left << std::setw(2 * total_depth + 5) << "[Index]"     //
           << std::left << std::setw(total_align_size + 10) << "[Kernel]"  //
           << "[Time]"                                                     //
           << std::right << std::setw(11) << "[Percent]" << std::endl;
    }

    ss << std::setw(2 * current_layer + 1) << "#"                                           //
       << std::setfill('0') << std::setw(2) << kernel_id                                    //
       << std::setfill('.') << std::setw(2 * (total_depth - current_layer) + 2) << ": "     //
       << std::setfill(' ') << std::left << std::setw(total_align_size + 10) << getName();  //
    if (is_timing) {
        ss << std::fixed << std::setprecision(3) << exec_time << "s";            //
        ss << std::right << std::setw(10) << std::fixed << std::setprecision(2)  //
           << 100 * exec_time / total_time << "%"                                //
           << std::endl;
    } else {
        ss << std::fixed << std::setprecision(3) << "--:--";                     //
        ss << std::right << std::setw(10) << std::fixed << std::setprecision(2)  //
           << "--:--"                                                            //
           << std::endl;
    }
    for (auto pkernel : _child_kernels) {
        ss << pkernel->generateInformationString(total_time, current_layer + 1, total_depth, total_align_size);
    }
    return ss.str();
}

void Kernel::setInputParam_impl(std::shared_ptr<Param> PM){};

void Kernel::setInputDataSet_impl(std::shared_ptr<DataSet> DS){};

Status& Kernel::initializeKernel_impl(Status& stat) { return stat; }

Status& Kernel::executeKernel_impl(Status& stat) { return stat; }

Status& Kernel::finalizeKernel_impl(Status& stat) { return stat; }

void Kernel::connectRelatedKernels(std::shared_ptr<Kernel>& ker) {
    bool found_kernel = false;
    for (auto&& pkernel : _all_kernels) {
        if (*ker == *pkernel) {
            found_kernel = true;
            break;
        }
    }
    if (!found_kernel) {
        _all_kernels.push_back(ker);
        for (auto&& pkernel : ker->_child_kernels) connectRelatedKernels(pkernel);
    }
}

std::map<std::string, Kernel*>& getDictOfKernels() {
    static std::map<std::string, Kernel*> static_dict_of_kernels;
    return static_dict_of_kernels;
}

std::vector<Kernel*>& Kernel::getKernels() {
    static std::vector<Kernel*> static_all_kernels;
    return static_all_kernels;
}
};  // namespace PROJECT_NS

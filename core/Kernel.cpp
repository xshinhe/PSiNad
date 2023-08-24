#include "Kernel.h"

namespace PROJECT_NS {

/**
 * @brief constructors
 */
Kernel::Kernel(const std::string& iname) : customized_name{iname} {
    TOTAL += 1;
    TOTAL_ROOT += 1;
    kernel_id      = TOTAL;
    max_align_size = name().size();
};

/**
 * @brief deconstructor.
 */
Kernel::~Kernel() {
    if (is_root) destory();
    TOTAL_ROOT -= 1;
}

/**
 * @brief build tree structure of the kernel
 */
Kernel& Kernel::push(std::shared_ptr<Kernel> ker, const bool& take_ownership) {
    _kernel_vector.push_back(ker);
    if (take_ownership) {
        _kernel_vector.back()->is_root = false;  // take the ownership
        TOTAL_ROOT -= 1;
    }

    depth = std::max(depth, ker->depth + 1);
    if (ker->name().size() > max_align_size) { max_align_size = ker->name().size(); }

    return *this;
}


/**
 * @brief  print information of the kernel
 */
std::string Kernel::scheme(double total_time, int current_layer, int total_depth, int total_align_size) {
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

    ss << std::setw(2 * current_layer + 1) << "#"                                                                   //
       << std::setfill('0') << std::setw(2) << kernel_id                                                            //
       << std::setfill('.') << std::setw(2 * (total_depth - current_layer) + 2) << ": "                             //
       << std::setfill(' ') << std::left << std::setw(total_align_size + 10) << name()                              //
       << std::fixed << std::setprecision(3) << exec_time << "s"                                                    //
       << std::right << std::setw(10) << std::fixed << std::setprecision(2) << 100 * exec_time / total_time << "%"  //
       << std::endl;
    for (auto pkernel : _kernel_vector) {
        ss << pkernel->scheme(total_time, current_layer + 1, total_depth, total_align_size);
    }
    return ss.str();
}

void Kernel::read_param_impl(Param* P){};  // overwritable param initializer

void Kernel::init_data_impl(DataSet* DS){};  // overwritable data reference

void Kernel::init_calc_impl(int stat){};  // overwritable data reference

int Kernel::exec_kernel_impl(int stat) { return 0; }  // overwritable run kernel

void Kernel::destory() {
    // automatically enabled by shared_ptr
    /**
    for (auto pkernel : _kernel_vector) {  // loop sub-kernels
        if (pkernel == nullptr || pkernel->is_root == true) continue;
        pkernel->destory();
        Kernel* ker = pkernel;
        delete ker;
        // pkernel = nullptr;
    }
     */
}

int Kernel::TOTAL      = 0;
int Kernel::TOTAL_ROOT = 0;
bool Kernel::BREAK     = false;

};  // namespace PROJECT_NS

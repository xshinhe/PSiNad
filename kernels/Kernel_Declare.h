#ifndef Kernel_Declare_H
#define Kernel_Declare_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * Dimension namespace locates read parameters for system size
 */
namespace Dimension {
extern int M;  ///< No. of MonteCarlo
extern int P;  ///< No. of Parallel
extern int N;  ///< No. of Nuclear phase space pairs
extern int F;  ///< No. of Fock State (electronic state)
extern int MP;
extern int PP;
extern int PN;
extern int PNN;
extern int PF;
extern int PFF;
extern int PNFF;
extern int NF;
extern int NN;
extern int FF;
extern int NFF;
extern int Fadd1;
// extern int nstep;
// extern int nsamp;
};  // namespace Dimension

/**
 * This class specifies the kernels that should establish a primary connection with dataset object.
 */
class Kernel_Declare final : public Kernel {
   public:
    Kernel_Declare(std::vector<std::shared_ptr<Kernel>> kers) : Kernel() {
        for (auto& ker : kers) { _ref_kernels.push_back(ker); }
    }

    inline virtual const std::string name() {
        std::stringstream ss;
        ss << "Kernel_Declare";
        for (auto& ker : _ref_kernels) ss << " #" << std::setfill('0') << std::setw(2) << ker->id();
        return ss.str();
    }

   private:
    std::vector<std::shared_ptr<Kernel>> _ref_kernels;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* DS);
};

/**
 * This class specifies the kernels that should be initialize at first.
 */
class Kernel_Initialize final : public Kernel {
   public:
    Kernel_Initialize(std::vector<std::shared_ptr<Kernel>> kers) : Kernel() {
        for (auto& ker : kers) { _ref_kernels.push_back(ker); }
    }

    inline virtual const std::string name() {
        std::stringstream ss;
        ss << "Kernel_Initialize";
        for (auto& ker : _ref_kernels) ss << " #" << std::setfill('0') << std::setw(2) << ker->id();
        return ss.str();
    }

   private:
    std::vector<std::shared_ptr<Kernel>> _ref_kernels;

    virtual void init_calc_impl(int stat = -1);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Declare_H
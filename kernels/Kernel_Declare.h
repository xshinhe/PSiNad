#ifndef Kernel_Declare_H
#define Kernel_Declare_H

#include "../core/Kernel.h"

#define TEST_V(NAME) std::cout << "test : " << #NAME << " = " << NAME << "\n";

namespace PROJECT_NS {

namespace Dimension {
extern int M;  // N_MonteCarlo;
extern int P;  // N_MonteCarlo;
extern int N;  // N_Nucl;
extern int F;  // N_Elec;
extern int MP;
extern int PN;
extern int NF;
extern int NN;
extern int FF;
extern int NFF;
extern int Fadd1;
// extern int nstep;
// extern int nsamp;
};  // namespace Dimension

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

};  // namespace PROJECT_NS

#endif  // Kernel_Declare_H
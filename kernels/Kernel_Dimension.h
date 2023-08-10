#ifndef Kernel_Dimension_H
#define Kernel_Dimension_H

#include "../core/Kernel.h"

#define TEST_V(NAME) std::cout << "test : " << #NAME << " = " << NAME << "\n";

namespace PROJECT_NS {
/**
 * @brief Kernel_Random manipulation of random engine and numbers
 */
class Kernel_Dimension final : public Kernel {
   public:
    virtual const std::string name() { return "Kernel_Dimension"; }

    static int M;  // N_MonteCarlo;
    static int P;  // N_MonteCarlo;
    static int N;  // N_Nucl;
    static int F;  // N_Elec;

    static int MP;
    static int PN;
    static int NF;
    static int NN;
    static int FF;
    static int NFF;
    static int Fadd1;

   private:
    virtual void read_param_impl(Param* PM);
};

};  // namespace PROJECT_NS

#endif  // Kernel_Dimension_H
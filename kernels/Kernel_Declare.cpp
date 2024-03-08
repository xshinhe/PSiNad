#include "Kernel_Declare.h"

namespace PROJECT_NS {

namespace Dimension {
int M     = 1;
int P     = 1;
int N     = 1;
int F     = 1;
int MP    = 1;
int PP    = 1;
int PN    = 1;
int PNN   = 1;
int PF    = 1;
int PFF   = 1;
int PNFF  = 1;
int NF    = 1;
int NN    = 1;
int FF    = 1;
int NFF   = 1;
int Fadd1 = 1;
};  // namespace Dimension

void Kernel_Declare::read_param_impl(Param* PM) {
    Dimension::M = PM->get<int>("M", LOC(), Dimension::M);
    Dimension::P = PM->get<int>("P", LOC(), Dimension::P);
    Dimension::N = PM->get<int>("N", LOC(), Dimension::N);
    Dimension::F = PM->get<int>("F", LOC(), Dimension::F);
    //
    // auxiliary dimension
    Dimension::FF    = Dimension::F * Dimension::F;
    Dimension::NFF   = Dimension::N * Dimension::FF;
    Dimension::NF    = Dimension::N * Dimension::F;
    Dimension::NN    = Dimension::N * Dimension::N;
    Dimension::PP    = Dimension::P * Dimension::P;
    Dimension::PN    = Dimension::P * Dimension::N;
    Dimension::PNN   = Dimension::P * Dimension::NN;
    Dimension::PF    = Dimension::P * Dimension::F;
    Dimension::PFF   = Dimension::P * Dimension::FF;
    Dimension::PNFF  = Dimension::P * Dimension::NFF;
    Dimension::MP    = Dimension::M * Dimension::P;
    Dimension::Fadd1 = Dimension::F + 1;

    for (auto& ker : _ref_kernels) ker->read_param(PM);
}

void Kernel_Declare::init_data_impl(DataSet* DS) {
    for (auto& ker : _ref_kernels) ker->init_data(DS);
}

void Kernel_Initialize::init_calc_impl(int stat) {
    for (auto& ker : _ref_kernels) ker->init_calc(stat);
}

};  // namespace PROJECT_NS
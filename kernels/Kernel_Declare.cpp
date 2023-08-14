#include "Kernel_Declare.h"

namespace PROJECT_NS {

namespace Dimension {
int M     = 1;  // N_MonteCarlo;
int P     = 1;  // N_MonteCarlo;
int N     = 1;  // N_Nucl;
int F     = 1;  // N_Elec;
int MP    = 1;
int PN    = 1;
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
    Dimension::PN    = Dimension::P * Dimension::N;
    Dimension::MP    = Dimension::M * Dimension::P;
    Dimension::Fadd1 = Dimension::F + 1;

    for (auto& ker : _ref_kernels) ker->read_param(PM);
}

void Kernel_Declare::init_data_impl(DataSet* DS) {
    for (auto& ker : _ref_kernels) ker->init_data(DS);
}

};  // namespace PROJECT_NS
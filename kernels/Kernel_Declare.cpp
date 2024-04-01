#include "Kernel_Declare.h"

#include "../core/vars_list.h"

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

    DATA::_M = PM->get<int>("M", LOC(), 1);
    DATA::_P = PM->get<int>("P", LOC(), 1);
    DATA::_N = PM->get<int>("N", LOC(), DATA::_N);
    DATA::_F = PM->get<int>("F", LOC(), DATA::_F);

    DATA::shape_M.static_build();
    DATA::shape_P.static_build();
    DATA::shape_N.static_build();
    DATA::shape_F.static_build();
    DATA::shape_Fadd1.static_build();

    DATA::shape_MP.static_build();
    DATA::shape_PP.static_build();
    DATA::shape_PN.static_build();
    DATA::shape_PNN.static_build();
    DATA::shape_PF.static_build();
    DATA::shape_PFF.static_build();
    DATA::shape_PNFF.static_build();
    DATA::shape_NF.static_build();
    DATA::shape_NN.static_build();
    DATA::shape_FF.static_build();
    DATA::shape_NFF.static_build();
    DATA::shape_NNFF.static_build();

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
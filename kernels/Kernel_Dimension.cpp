#include "Kernel_Dimension.h"

namespace PROJECT_NS {

int Kernel_Dimension::M = 1;  // N_MonteCarlo;
int Kernel_Dimension::P = 1;  // N_MonteCarlo;
int Kernel_Dimension::N = 1;  // N_Nucl;
int Kernel_Dimension::F = 1;  // N_Elec;

int Kernel_Dimension::MP    = 1;
int Kernel_Dimension::PN    = 1;
int Kernel_Dimension::NF    = 1;
int Kernel_Dimension::NN    = 1;
int Kernel_Dimension::FF    = 1;
int Kernel_Dimension::NFF   = 1;
int Kernel_Dimension::Fadd1 = 1;

void Kernel_Dimension::read_param_impl(Param* PM) {
    M = PM->get<int>("M", LOC(), M);
    P = PM->get<int>("P", LOC(), P);
    N = PM->get<int>("N", LOC(), N);
    F = PM->get<int>("F", LOC(), F);

    // first read "model"
    auto& j0 = (*(PM->pjson()));
    if (j0.count("model_param") > 0) {
        auto& j1 = j0["model_param"];
        if (j1.count("M") > 0) M = j1["M"].as_integer();
        if (j1.count("P") > 0) M = j1["P"].as_integer();
        if (j1.count("N") > 0) M = j1["N"].as_integer();
        if (j1.count("F") > 0) M = j1["F"].as_integer();
    }

    // auxiliary dimension
    FF    = F * F;
    NFF   = N * FF;
    NF    = N * F;
    NN    = N * N;
    PN    = P * N;
    MP    = M * P;
    Fadd1 = F + 1;
}

};  // namespace PROJECT_NS
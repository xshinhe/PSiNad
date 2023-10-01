#include "ns.h"

#include <iostream>

// #include "types.h"

#define DEF_VARIABLE(TYPE, name, doc) VARIABLE<TYPE> name = VARIABLE<TYPE>(#name, doc);

namespace PROJECT_NS {


namespace DATASET {

DEF_VARIABLE(int, dimension::M, "Dimension of Monte Carlo Size");
DEF_VARIABLE(int, dimension::P, "Parallel");
DEF_VARIABLE(double, dimension::x, "coordinate");
DEF_VARIABLE(double, integrator::baoab::x, "coordinate 2");

};  // namespace DATASET

};  // namespace PROJECT_NS

int main() {
    std::cout << PROJECT_NS::DATASET::dimension::M.help() << "\n";
    std::cout << PROJECT_NS::DATASET::dimension::M.name() << "\n";

    for (auto& i : PROJECT_NS::DATASET::_VAR_LIST) {
        std::cout << ((PROJECT_NS::DATASET::VARIABLE<int>*) i)->name() << " : "
                  << ((PROJECT_NS::DATASET::VARIABLE<int>*) i)->help() << "\n";
    }

    return 0;
}
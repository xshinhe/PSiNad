#include "ns.h"

#include <iostream>

#include "DataSet_test.h"

// #include "types.h"

#define DEF_VARIABLE(TYPE, name, shape, doc) VARIABLE<TYPE> name = VARIABLE<TYPE>(#name, shape, doc);
#define OPENDF_DECLARE_VARIABLE(type, name) \
    namespace name {                        \
    VARIABLE<type> var;                     \
    };

#define OPENDF_DEFINE_VARIABLE(type, name, shape, doc) \
    namespace name {                                   \
    VARIABLE<type> var(#name, shape, doc);             \
    };

namespace PROJECT_NS {

namespace apple::banana {
int c = 0;
};  // namespace apple::banana

namespace DATASET {
std::vector<void*> _VAR_LIST;

namespace DIM {
DSize P;
DSize N;
DSize F;
DSize M;
Shape PN({&P, &N});
};  // namespace DIM

DEF_VARIABLE(double, dimension::x, &DIM::PN, "coordinate");
DEF_VARIABLE(double, integrator::baoab::x, &DIM::PN, "coordinate 2");


OPENDF_DEFINE_VARIABLE(double, integrator::x, &DIM::PN, "coordinates of integrator");
OPENDF_DEFINE_VARIABLE(double, integrator::p, &DIM::PN, "coordinates of integrator");


};  // namespace DATASET

};  // namespace PROJECT_NS

using namespace PROJECT_NS;

int main() {
    DataSet DS;
    for (auto& i : PROJECT_NS::DATASET::_VAR_LIST) {
        std::cout << ((PROJECT_NS::DATASET::VARIABLE<int>*) i)->name() << " : "
                  << ((PROJECT_NS::DATASET::VARIABLE<int>*) i)->help() << "\n";
    }
    std::cout << DATASET::integrator::x::var.name() << "\n";

    double* x = DS.reg(DATASET::integrator::x::var);



    return 0;
}

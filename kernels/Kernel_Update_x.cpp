#include "Kernel_Update_x.h"

#include "Kernel_Declare.h"

#define ARRAY_SHOW(_A, _n1, _n2)                                                            \
    ({                                                                                      \
        std::cout << #_A << " = np.array([\n";                                              \
        int _idxA = 0;                                                                      \
        for (int _i = 0; _i < (_n1); ++_i) {                                                \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(8) << (_A)[_idxA++] << ","; \
            std::cout << std::endl;                                                         \
        }                                                                                   \
        { std::cout << "])\n"; }                                                            \
    })

namespace PROJECT_NS {

void Kernel_Update_x::init_data_impl(DataSet* DS) {
    dt_ptr          = DS->def<kids_real>("iter.dt");
    x               = DS->def<kids_real>("integrator.x", Dimension::PN);
    p               = DS->def<kids_real>("integrator.p", Dimension::PN);
    m               = DS->def<kids_real>("integrator.m", Dimension::PN);
    minv            = DS->def<kids_real>("integrator.minv", Dimension::PN);
    kids_real* mass = DS->def<kids_real>("model.mass", Dimension::N);
    frez_ptr        = DS->def<kids_bool>("iter.frez");
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* m    = this->m + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        for (int j = 0; j < Dimension::N; ++j) {
            m[j]    = mass[j];
            minv[j] = 1 / m[j];
        }
    }
}

int Kernel_Update_x::exec_kernel_impl(int stat) {
    if (frez_ptr[0]) return 0;
    for (int i = 0; i < Dimension::PN; ++i) x[i] += p[i] * minv[i] * scale * dt_ptr[0];
    return 0;
}

};  // namespace PROJECT_NS

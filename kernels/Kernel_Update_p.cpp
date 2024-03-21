#include "Kernel_Update_p.h"

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

void Kernel_Update_p::init_data_impl(DataSet* DS) {
    dt_ptr   = DS->def<kids_real>("iter.dt");
    f        = DS->def<kids_real>("integrator.f", Dimension::PN);
    p        = DS->def<kids_real>("integrator.p", Dimension::PN);
    minv     = DS->def<kids_real>("integrator.minv", Dimension::PN);
    Ekin     = DS->def<kids_real>("integrator.Ekin", Dimension::P);
    frez_ptr = DS->def<bool>("iter.frez");
}

void Kernel_Update_p::init_calc_impl(int stat) {
    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;
        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
    }
}

int Kernel_Update_p::exec_kernel_impl(int stat) {
    if (frez_ptr[0]) return 0;

    for (int iP = 0; iP < Dimension::P; ++iP) {
        kids_real* f    = this->f + iP * Dimension::N;
        kids_real* p    = this->p + iP * Dimension::N;
        kids_real* minv = this->minv + iP * Dimension::N;
        kids_real* Ekin = this->Ekin + iP;

        //////////////////////////////////////////////
        Ekin[0] = 0.0e0;
        for (int i = 0; i < Dimension::N; ++i) {
            p[i] -= f[i] * scale * dt_ptr[0];
            Ekin[0] += 0.5e0 * p[i] * p[i] * minv[i];
        }
    }
    return 0;
}

};  // namespace PROJECT_NS
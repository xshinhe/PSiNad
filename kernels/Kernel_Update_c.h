#ifndef Kernel_Update_c_H
#define Kernel_Update_c_H

#include "../core/Kernel.h"

namespace kids {

class Kernel_Update_c final : public Kernel {
   public:
    Kernel_Update_c(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_c"; }

   private:
    kids_complex* Udt;  ///< short time propagator
    kids_complex* U;    ///< full propagator along classical path approximation (CPA)

    ///< solve Diabatic propagator
    kids_real* E;  ///< Eigenvalue for diabatic V
    kids_real* T;  ///< Eigenvector for diabatic V

    ///< solve Adiabatic propagator
    kids_real* L;     ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
    kids_complex* R;  ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

    kids_complex* invexpidiagdt;  ///< temporary variables

    double scale, *dt_ptr;

    virtual void init_data_impl(DataSet* S);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace kids


#endif  // Kernel_Update_c_H

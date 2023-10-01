#ifndef Kernel_Update_c_H
#define Kernel_Update_c_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_c final : public Kernel {
   public:
    Kernel_Update_c(double scale) : Kernel(), scale{scale} {};

    inline virtual const std::string name() { return "Kernel_Update_c"; }

   private:
    num_complex* Udt;  ///< short time propagator
    num_complex* U;    ///< full propagator along classical path approximation (CPA)

    ///< solve Diabatic propagator
    num_real* E;  ///< Eigenvalue for diabatic V
    num_real* T;  ///< Eigenvector for diabatic V

    ///< solve Adiabatic propagator
    num_real* L;     ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
    num_complex* R;  ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

    num_complex* invexpidiagdt;  ///< temporary variables

    double scale;
    double dt, sdt;

    virtual void read_param_impl(Param* PM);

    virtual void init_data_impl(DataSet* S);

    virtual int exec_kernel_impl(int stat = -1);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_c_H

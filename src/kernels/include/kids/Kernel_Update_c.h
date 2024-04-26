#ifndef Kernel_Update_c_H
#define Kernel_Update_c_H

#include "kids/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_c final : public Kernel {
   public:
    Kernel_Update_c(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

   private:
    kids_complex* Udt;  ///< short time propagator
    kids_complex* U;    ///< full propagator along classical path approximation (CPA)

    ///< solve Diabatic propagator
    kids_real* E;  ///< Eigenvalue for diabatic V
    kids_real* T;  ///< Eigenvector for diabatic V

    ///< solve Adiabatic propagator
    kids_real*    L;  ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
    kids_complex* R;  ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

    kids_complex* invexpidiagdt;  ///< temporary variables

    double     scale, *dt_ptr;
    kids_bint* succ_ptr;
    kids_bint* frez_ptr;

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    virtual Status& executeKernel_impl(Status& stat);
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_c_H

/**@file        Kernel_Update_U.h
 * @brief       this file provides Kernel_Update_U class for update U propagator
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is a product of academic research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @par revision:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-02  <td> Initial version. Added detailed commentary by ChatGPT.
 * </table>
 **********************************************************************************
 */

#ifndef Kernel_Update_U_H
#define Kernel_Update_U_H

#include "psnd/Kernel.h"

namespace PROJECT_NS {

class Kernel_Update_U final : public Kernel {
   public:
    Kernel_Update_U(double scale) : Kernel(), scale{scale} {};

    virtual const std::string getName();

    virtual int getType() const;

   private:
    bool only_adjust;
    bool enable_update_c;
    bool enable_update_cset;
    bool enable_update_rho_ele;
    bool enable_update_rho_nuc;
    bool enable_update_rho_dual;

    span<psnd_complex> Udt;  ///< short time propagator
    span<psnd_complex> U;    ///< full propagator along classical path approximation (CPA)

    ///< solve Diabatic propagator
    span<psnd_real> eig, dE;    ///< Eigenvalue for diabatic V
    span<psnd_real> T, T_init;  ///< Eigenvector for diabatic V

    ///< solve Adiabatic propagator
    span<psnd_real>    lam;  ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
    span<psnd_complex> R;    ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

    span<psnd_complex> invexpidiagdt;  ///< temporary variables

    span<psnd_complex> c, c_init;
    span<psnd_complex> cset, cset_init;
    span<psnd_complex> rho_ele, rho_ele_init;
    span<psnd_complex> rho_nuc, rho_nuc_init;
    span<psnd_complex> rho_dual, rho_dual_init;

    span<psnd_real>    mono, monodt;
    span<psnd_complex> MFFtmp1, MFFtmp2;

    psnd_real       scale;
    span<psnd_real> dt;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status &initializeKernel_impl(Status &stat);

    virtual Status &executeKernel_impl(Status &stat);

    void update_monodromy();
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_U_H

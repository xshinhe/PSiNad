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
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
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

#include "kids/Kernel.h"

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

    span<kids_complex> Udt;  ///< short time propagator
    span<kids_complex> U;    ///< full propagator along classical path approximation (CPA)

    ///< solve Diabatic propagator
    span<kids_real> eig, dE;    ///< Eigenvalue for diabatic V
    span<kids_real> T, T_init;  ///< Eigenvector for diabatic V

    ///< solve Adiabatic propagator
    span<kids_real>    lam;  ///< Eigenvalue for adiabatic effective Hamiltonian Heff = Eδ - id*P/M
    span<kids_complex> R;    ///< Eigenvector for adiabatic effective Hamiltonian Heff = Eδ - id*P/M

    span<kids_complex> invexpidiagdt;  ///< temporary variables

    span<kids_complex> c, c_init;
    span<kids_complex> cset, cset_init;
    span<kids_complex> rho_ele, rho_ele_init;
    span<kids_complex> rho_nuc, rho_nuc_init;
    span<kids_complex> rho_dual, rho_dual_init;

    span<kids_real>    mono, monodt;
    span<kids_complex> MFFtmp1, MFFtmp2;

    kids_real       scale;
    span<kids_real> dt;

    virtual void setInputParam_impl(std::shared_ptr<Param> PM);

    virtual void setInputDataSet_impl(std::shared_ptr<DataSet> DS);

    virtual Status &initializeKernel_impl(Status &stat);

    virtual Status &executeKernel_impl(Status &stat);

    void update_monodromy();
};


};  // namespace PROJECT_NS


#endif  // Kernel_Update_U_H

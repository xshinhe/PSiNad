#ifndef BATH_H
#define BATH_H
#include <map>

#include "../core/Kernel.h"

namespace PROJECT_NS {

namespace BathPolicy {

enum _enum_coup { SpinBosonCoup, SiteExcitonCoup, GeneralCoup };
const std::map<std::string, int> _dict_coup = {
    {"#sb", SpinBosonCoup},    // sigma_z is coulped
    {"#se", SiteExcitonCoup},  // {|n><n|,n=1,...k; k<=F} is coupled (k can be less than F)
    {"###", GeneralCoup}       // coupling matrix is read from "coup.dat" file
};

enum _enum_lamb { lambUsed, alphaUsed, etaUsed, ergUsed };
const std::map<std::string, int> _dict_lamb = {
    {"lambda", lambUsed},  // lambda: characherized in Debye spectrum
    {"alpha", alphaUsed},  // alpha: kondo parameter in Ohmic spectrum
    {"eta", etaUsed},      // eta: friction parameter in Brownian spectrum
    {"erg", ergUsed}       // erg: renormalization energy of a spectrum
};

enum _enum_spec { DebyeSpec, OhmicSpec, ClosureSpec, HuangRhysSpec, GFactorSpec, ReadSpec };
const std::map<std::string, int> _dict_spec = {
    {"debye", DebyeSpec},          // Debye spectrum (continuous spectrum)
    {"ohmic", OhmicSpec},          // Ohmic spectrum (continuous spectrum)
    {"closure", ClosureSpec},      // Any closure decomposition of TCF of the bath, it will read from "spec.dat"
    {"#b850", HuangRhysSpec},      // intrinsic B850 spectrum (by Huang-Rhys format), see hamilontian_data.h
    {"#pbi", HuangRhysSpec},       // intrinsic PBI spectrum (by Huang-Rhys format), see hamilontian_data.h
    {"huangrhys", HuangRhysSpec},  // general Huang-Rhys format (Discretized spectrum)
    {"#rub", GFactorSpec},         // intrinsic rubenzene spectrum (by g-factor format), see hamilontian_data.h
    {"gfactor", GFactorSpec},      // general g-factor format (Discretized spectrum)
    {"read", ReadSpec}             // read format (omegas, coeffs)
};
};  // namespace BathPolicy

class Bath final : public Kernel {
   public:
    int Discretization(double* W_arr, double* C_arr, const int& Nb);

    static inline double J_Ohmic(const double& w);

    static inline double J_Debye(const double& w);

    static inline double J_Refit(const double& w, double* W_arr, double* C_arr, const int& Nb);

    static inline double J(const double& w, double* W_arr = nullptr, double* C_arr = nullptr, const int& Nb = 0);

    static inline double fun_Cw_Re(const double& w, double* W_arr, double* C_arr, const int& Nb);

    static inline int fun_Cw(num_complex* Cw_arr, double* w, const int& Nw, double* W_arr, double* C_arr,
                             const int& Nb);

   private:
    Option lambda_option;
    Option coupling_option;
    Option spectrum_option;

    // description for each bath
    num_real omegac;  ///< characheristic frequency
    num_real lambda;  ///< (proportional to) reorganization energy
    num_real beta;    ///<  1/(k_B * Temperature)

    int nbath, F;
    num_real* Q;

    virtual void read_param_impl(Param* P) {
        // read omegac
        omegac = P->get<double>("omegac", LOC(), phys::energy_d, 1.0f);
        CHECK_GT(omegac, 0);

        lambda_option.flag = P->get<std::string>("lambda_flag", LOC(), phys::none_d, "lambda");
        lambda_option.type = lambda_policy::_dict.at(lambda_option.flag);

        switch (lambda_option.type) {
            case lambda_policy::lambUsed: {
                double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
                lambda          = strength;
                break;
            }
            case lambda_policy::alphaUsed: {
                double strength = P->get<double>("strength", LOC(), phys::none_d, 1.0f);  // dimensionless
                lambda          = 0.5f * omegac * strength;
                break;
            }
            case lambda_policy::etaUsed: {
                double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
                lambda          = 0.5f * strength;
                break;
            }
            case lambda_policy::ergUsed: {
                double strength = P->get<double>("strength", LOC(), phys::energy_d, 1.0f);
                lambda          = 0.25f * strength;
                break;
            }
            default:
                LOG(FATAL) << "Unknown st_flag";
        }
        CHECK_GT(lambda, 0);

        // read temperature
        double temp = P->get<double>("temp", LOC(), phys::temperature_d, 1.0f);
        beta        = 1.0f / (phys::au::k * temp);  // never ignore k_Boltzman
        CHECK_GT(beta, 0);

        // read coupling flag
        coupling_option.flag = P->get<std::string>("coupling_flag", LOC(), phys::none_d, "#se");
        coupling_option.type = coupling_policy::_dict.at(coupling_option.flag);

        // read spectrum flag
        spectrum_option.flag = P->get<std::string>("spectrum_flag", LOC(), phys::none_d, "read");
        spectrum_option.type = spectrum_policy::_dict.at(spectrum_option.flag);


        // allocate Q tensor (the coupling ternsor) & initialization
        nbath = nbath_in, F = F_in;
        ALLOCATE_PTR_TO_VECTOR(Q, nbath * F * F);
        memset(Q, 0, nbath * F * F * sizeof(num_real));

        switch (coup_type) {
            case BathPolicy::SpinBosonCoup: {
                CHECK_EQ(nbath, 1);
                CHECK_EQ(F, 2);                                       // nbath must be 1, F must be 2
                Q[0] = 1.0f, Q[1] = 0.0f, Q[2] = 0.0f, Q[3] = -1.0f;  // sigmaz
                break;
            }
            case BathPolicy::SiteExcitonCoup: {
                // always be careful with nbath and F
                CHECK_LE(nbath, F);  // nbath should equal or less than F. (a ground state can be placed in last)
                for (int i = 0, idx = 0; i < nbath; ++i) {
                    for (int j = 0; j < F; ++j) {
                        for (int k = 0; k < F; ++k, ++idx) Q[idx] = (i == j && i == k) ? 1.0f : 0.0f;
                    }
                }
                break;
            }
            case BathPolicy::GeneralCoup:
            default: {
                std::ifstream ifs(coup_flag);
                num_real tmp;
                for (int i = 0; i < nbath * F * F; ++i)
                    if (ifs >> tmp) Q[i] = tmp;
                ifs.close();
            }
        }
    }

    virtual void init_data_impl(State* S) {
        // external
        x = S->reg<double>("integrator.x", size);
        f = S->reg<double>("integrator.f", size);
        // shared
        w = S->reg<double>("model.w", size);
        m = S->reg<double>("model.m", size);
        for (int i = 0; i < size; ++i) m[i] = 1, w[i] = 1;
    }

    virtual int exec_kernel_impl(int stat = -1) {
        for (int i = 0; i < size; ++i) f[i] = m[i] * w[i] * w[i] * x[i];
        return 0;
    }
};
};  // namespace PROJECT_NS

#endif  // BATH_H

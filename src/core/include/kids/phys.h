#ifndef PHYS_H
#define PHYS_H

#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>

// iso standard:
// C++98   #define __cplusplus 199711
// C++11   #define __cplusplus 201103
// C++14   #define __cplusplus 201402
// C++17   #define __cplusplus 201703
// C++20   #define __cplusplus 202002

#if (__cplusplus < 201402L)
/** make this header is compatible to c++11 standard */
#define CONSTTYPE const
#define CONSTEXPR_DECOR
#else
#define FULLY_COMPILE_TIME_OPT
#define CONSTTYPE constexpr
#define CONSTEXPR_DECOR constexpr
#endif  // (__cplusplus < 201402L)

// using real_precision = kids_real;
using real_precision = double;

namespace phys {  // phys::

namespace math {  // phys::math::

constexpr real_precision eu     = 2.718281828459045235360287;   //!< Euler'Constant
constexpr real_precision pi     = 3.141592653589793238462643L;  //!< pi
constexpr real_precision twopi  = 6.283185307179586476925287L;
constexpr real_precision halfpi = 1.570796326794896619231321L;
constexpr real_precision eps8   = 1.0E-8L;
constexpr real_precision eps16  = 1.0E-16L;
constexpr real_precision eps32  = 1.0E-32L;

constexpr real_precision sqrttwo  = 1.414213562373095048801689L;
constexpr real_precision sqrthalf = 0.707106781186547524400844L;

constexpr std::complex<real_precision> im(0.0L, 1.0L);  //!< Imaginary Unit
constexpr std::complex<real_precision> iu(1.0L, 0.0L);
constexpr std::complex<real_precision> iz(0.0L, 0.0L);

};  // namespace math

struct unit_error : public std::runtime_error {
    unit_error(std::string const text) : std::runtime_error(text) {}
};

////////////////////////////////////////////////////////////////////////////////

#if (__cplusplus < 201703L)
/**
 * @brief compile-time array operations realized in c++11 standard
 */

template <int...>
struct seq {};
template <int N, int... S>
struct gens : gens<N - 1, N - 1, S...> {};
template <int... S>
struct gens<0, S...> {
    using type = seq<S...>;
};
template <typename T, std::size_t N>
constexpr T array_add_n(const std::array<T, N> a, const std::array<T, N> b, int n) {
    return a[n] + b[n];
}
template <typename T, std::size_t N, int... S>
constexpr std::array<T, N> array_add_impl(const std::array<T, N> a, const std::array<T, N> b, seq<S...>) {
    return std::array<T, N>{{array_add_n(a, b, S)...}};
}
template <typename T, std::size_t N>
constexpr std::array<T, N> array_add(const std::array<T, N> a, const std::array<T, N> b) {
    return array_add_impl(a, b, typename gens<N>::type());
}

template <typename T, std::size_t N>
constexpr T array_minus_n(const std::array<T, N> a, const std::array<T, N> b, int n) {
    return a[n] - b[n];
}
template <typename T, std::size_t N, int... S>
constexpr std::array<T, N> array_minus_impl(const std::array<T, N> a, const std::array<T, N> b, seq<S...>) {
    return std::array<T, N>{{array_minus_n(a, b, S)...}};
}
template <typename T, std::size_t N>
constexpr std::array<T, N> array_minus(const std::array<T, N> a, const std::array<T, N> b) {
    return array_minus_impl(a, b, typename gens<N>::type());
}

template <typename T, std::size_t N>
constexpr T array_scale_n(const std::array<T, N> a, T b_val, int n) {
    return a[n] * b_val;
}
template <typename T, std::size_t N, int... S>
constexpr std::array<T, N> array_scale_impl(const std::array<T, N> a, const T b, seq<S...>) {
    return std::array<T, N>{{array_scale_n(a, b, S)...}};
}
template <typename T, std::size_t N>
constexpr std::array<T, N> array_scale(const std::array<T, N> a, T b) {
    return array_scale_impl(a, b, typename gens<N>::type());
}

#else  // (if c++17 is support it will be much easier)

template <typename T, std::size_t N>
constexpr std::array<T, N> array_add(const std::array<T, N> a, const std::array<T, N> b) {
    std::array<T, N> res{};
    for (int i = 0; i < N; ++i) res[i] = a[i] + b[i];
    return res;
}
template <typename T, std::size_t N>
constexpr std::array<T, N> array_minus(const std::array<T, N> a, const std::array<T, N> b) {
    std::array<T, N> res{};
    for (int i = 0; i < N; ++i) res[i] = a[i] - b[i];
    return res;
}
template <typename T, std::size_t N>
constexpr std::array<T, N> array_scale(const std::array<T, N> a, T b) {
    std::array<T, N> res{};
    for (int i = 0; i < N; ++i) res[i] = a[i] * b;
    return res;
}
#endif  // (__cplusplus < 201703L)
////////////////////////////////////////////////////////////////////////////////

template <typename T, std::size_t N>
class dimensions {
   public:
    using type = T;
    std::array<T, N> _data;

    constexpr dimensions() : _data{0} {};

    constexpr dimensions(const std::array<T, N>& a) : _data{a} {};

    // for aggregate initialization, please using dimensions{{....}};

    inline constexpr bool operator==(const dimensions& rdim) const { return (_data == rdim._data); }
    inline constexpr bool operator!=(const dimensions& rdim) const { return !(_data == rdim._data); }
    inline constexpr bool operator<(const dimensions& rdim) const { return (_data < rdim._data); }
    inline constexpr bool operator>(const dimensions& rdim) const { return (_data > rdim._data); }

    inline constexpr dimensions operator/(const dimensions rdim) const {
        return dimensions{array_minus(_data, rdim._data)};
    }
    inline constexpr dimensions operator*(const dimensions rdim) const {
        return dimensions{array_add(_data, rdim._data)};
    }
    inline constexpr dimensions power(const T index) const { return dimensions{array_scale(_data, index)}; }
    inline std::string          to_string() const {
        std::ostringstream os;
        os << "<";
        for (auto i : _data) os << i << ", ";
        os << ">";
        return os.str();
    }
};

/**
 * @brief dimension7 is provided as compile-time dimensional tools
 */
const int                                           dimension7_size = 7;
typedef dimensions<real_precision, dimension7_size> dimension7;
enum dimension7_type {
    _L,  // length dimension7
    _T,  // time dimension7
    _M,  // mass dimension7
    _I,  // electric current dimension7 (_C is for charge as alternate)
    _Q,  // thermodynamic temperature dimension7 (Q actually means \Theta)
    _N,  // amount of substance dimension7
    _J,  // luminous intensity dimension7
};

inline constexpr real_precision reduce_l_nonzero(const dimension7 dim) {
    return (int{dim._data[0] != 0L} + int{dim._data[1] != 0L} + int{dim._data[2] != 0L} + int{dim._data[3] != 0L} +
            int{dim._data[4] != 0L} + int{dim._data[5] != 0L} + int{dim._data[6] != 0L});
}
inline constexpr real_precision reduce_l_energy(const dimension7 dim) {
    return -dim._data[_L] - dim._data[_T] + dim._data[_M] + dim._data[_Q];
}

/** \name base dimension7 */
/// @{
constexpr dimension7 dimensionless_d{{}};                                 ///< [1]
constexpr dimension7 length_d{{1, 0, 0, 0, 0, 0, 0}};                     ///< [L]
constexpr dimension7 time_d{{0, 1, 0, 0, 0, 0, 0}};                       ///< [T]
constexpr dimension7 mass_d{{0, 0, 1, 0, 0, 0, 0}};                       ///< [M]
constexpr dimension7 electric_current_d{{0, 0, 0, 1, 0, 0, 0}};           ///< [I]
constexpr dimension7 thermodynamic_temperature_d{{0, 0, 0, 0, 1, 0, 0}};  ///< [Q]
constexpr dimension7 amount_of_substance_d{{0, 0, 0, 0, 0, 1, 0}};        ///< [N]
constexpr dimension7 luminous_intensity_d{{0, 0, 0, 0, 0, 0, 1}};         ///< [J]

constexpr dimension7 current_d     = electric_current_d;
constexpr dimension7 temperature_d = thermodynamic_temperature_d;
constexpr dimension7 amount_d      = amount_of_substance_d;
constexpr dimension7 none_d        = dimensionless_d;
/// @}

/** \name derived (L) dimension7 */
/// @{
constexpr dimension7 distance_d   = length_d;
constexpr dimension7 wavelength_d = length_d;
constexpr dimension7 wave_number_d{{-1, 0, 0}};  ///< [L^-1]
constexpr dimension7 area_d{{2, 0, 0}};          ///< [L^2]
constexpr dimension7 volume_d{{3, 0, 0}};        ///< [L^3]
/// @}

/** \name derived (T) dimension7 */
/// @{
constexpr dimension7 frequency_d{{0, -1, 0}};             ///< [T^-1]
constexpr dimension7 angular_velocity_d{{0, -1, 0}};      ///< [T^-1]
constexpr dimension7 angular_acceleration_d{{0, -2, 0}};  ///< [T^-2]
constexpr dimension7 activity_of_a_nuclide_d = frequency_d;
/// @}

/** \name derived (L, T) dimension7 */
/// @{
constexpr dimension7 speed_d{{1, -1, 0}};                       ///< [L*T^-1]
constexpr dimension7 acceleration_d{{1, -2, 0}};                ///< [L*T^-2]
constexpr dimension7 jerk_d{{1, -3, 0}};                        ///< [L*T^-3]
constexpr dimension7 jounce_d{{1, -4, 0}};                      ///< [L*T^-4]
constexpr dimension7 crackle_d{{1, -5, 0}};                     ///< [L*T^-5]
constexpr dimension7 pop_d{{1, -6, 0}};                         ///< [L*T^-6]
constexpr dimension7 absement_d{{1, 1, 0}};                     ///< [L*T]
constexpr dimension7 area_flow_rate_d{{2, -1, 0}};              ///< [L^2*T^-1]
constexpr dimension7 volume_flow_rate_d{{3, -1, 0}};            ///< [L^3*T^-1]
constexpr dimension7 kinematic_viscosity_d = area_flow_rate_d;  ///< [L^2*T^-1] = viscosity / density
constexpr dimension7 thermal_diffusivity_d =
    area_flow_rate_d;  ///< [L^2*T^-1] = thermal_conductivity / (specific_heat_capacity * density)
constexpr dimension7 specific_energy_d{{2, -2, 0}};          ///< [L^2/T^2] (count) energy per mass
constexpr dimension7 dose_equivalent_d = specific_energy_d;  ///< [L^2/T^2] (radiation) energy per mass
constexpr dimension7 absorbed_dose_d   = specific_energy_d;  ///< [L^2/T^2] (radiation) energy per mass
constexpr dimension7 absorbed_dose_rate_d{{2, -3, 0}};       ///< [L^2/T^3] (radiation) power per mass
constexpr dimension7 substance_permeability_d{{-1, 1, 0}};   ///< [L^-1*T]
/// @}

/** \name derived (L, T, M) dimension7 */
/// @{
constexpr dimension7 inertia_d{{2, 0, 1}};                  ///< [M*L^2]
constexpr dimension7 mass_line_density_d{{-1, 0, 1}};       ///< [M/L] mass per line
constexpr dimension7 mass_area_density_d{{-2, 0, 1}};       ///< [M/L^2] mass per area
constexpr dimension7 mass_density_d{{-3, 0, 1}};            ///< [M/L^3] mass per volume
constexpr dimension7 specific_volume_d{{3, 0, -1}};         ///< [L^3/M] volume per mass
constexpr dimension7 mass_flow_rate_d{{0, -1, 1}};          ///< [M/T] mass per time
constexpr dimension7 mass_flow_acceleration_d{{0, -2, 1}};  ///< [M/T^2] mass per per time
constexpr dimension7 mass_flow_jerk_d{{0, -3, 1}};          ///< [M/T^3] mass per per per time

constexpr dimension7 force_d{{1, -2, 1}};           ///< [M*L/T^2] mass times acceleration
constexpr dimension7 momentum_d{{1, -1, 1}};        ///< [M*L/T] force integrate time
constexpr dimension7 energy_d{{2, -2, 1}};          ///< [M*L^2/T^2] force integrate length
constexpr dimension7 moment_of_force_d = energy_d;  ///< [M*L^2/T^2] force cross length
constexpr dimension7 torque_d          = moment_of_force_d;
constexpr dimension7 angular_momentum_d{{2, -1, 1}};  ///< [M*L^2/T] torque integrate time
constexpr dimension7 action_d = angular_momentum_d;   ///< [M*L^2/T] energy integrate time

constexpr dimension7 inv_ener_d{{-2, 2, -1}};  ///< [M^-1*L^-2*T^2], inversed energy, such as 1/(kB * T)

constexpr dimension7 power_d{{2, -3, 1}};                           ///< [M*L^2/T^3] energy per time
constexpr dimension7 energy_density_d{{-1, -2, 1}};                 ///< [M/L/T^2] energy per volume
constexpr dimension7 pressure_d        = energy_density_d;          ///< [M/L/T^2] energy per volume = force per area
constexpr dimension7 surface_tension_d = mass_flow_acceleration_d;  ///< [M/T^2] energy per area
constexpr dimension7 energy_line_density_d{{1, -2, 1}};             ///< [M/L/T^2] energy per line
constexpr dimension7 power_density_d{{-1, -3, 1}};                  ///< [M/L/T^3] power per volume
constexpr dimension7 power_area_density_d = mass_flow_jerk_d;       ///< [M/T^3] power per area

/**
 * @note dynamic_viscosity (mu) is differ from kinematic_viscosity kinematic_viscosity_d (nu):
 *
 *      kinematic_viscosity = dynamic_viscosity / density
 */
constexpr dimension7 dynamic_viscosity_d{{-1, -1, 1}};  ///< [M/L/T] force / (area * gradient(velocity))


constexpr dimension7 heat_flow_rate_d         = power_d;                   ///< [M/L/T^3] (heat) energy per time
constexpr dimension7 heat_density_d           = mass_flow_acceleration_d;  ///< [M/T^2] (heat flow) energy per area
constexpr dimension7 heat_density_flow_rate_d = power_area_density_d;  ///< [M/T^3] (heat flow) energy per area per time
constexpr dimension7 heat_flux_density_d      = power_area_density_d;  ///< [M/T^3] (heat) energy per time per area

constexpr dimension7 radiant_intensity_d = power_d;               ///< [M/L/T^3] (radiation) energy per time
constexpr dimension7 radiance_d          = power_area_density_d;  ///< [M/T^3] (radiation) power per area
constexpr dimension7 irradiance_d        = power_area_density_d;  ///< [M/T^3] (radiation) power per area
/// @}

/** \name derived (L, T, M, I) dimension7 */
/// @{
constexpr dimension7 current_density_d{{-2, 0, 0, 1}};               ///< [I/L^2] current per area
constexpr dimension7 electric_charge_d{{0, 1, 0, 1}};                ///< [I*T] current integrate time
constexpr dimension7 electric_charge_density_d{{-3, 1, 0, 1}};       ///< [I/L^3*T] charge per volume
constexpr dimension7 electric_area_charge_density_d{{-2, 1, 0, 1}};  ///< [I/L^2*T] charge per area
constexpr dimension7 electric_line_charge_density_d{{-1, 1, 0, 1}};  ///< [I/L*T] charge per line
constexpr dimension7 electric_dipole_moment_d{{1, 1, 0, 1}};         ///< [I*L*T] charge times length

constexpr dimension7 electric_flux_density_d       = electric_area_charge_density_d;  ///< [I/L^2*T]
constexpr dimension7 electric_displacement_field_d = electric_area_charge_density_d;  ///< [I/L^2*T], D
constexpr dimension7 electric_polarization_field_d =
    electric_area_charge_density_d;  ///< [I/L^2*T], P = dipole moment pe volume

constexpr dimension7 magnetic_moment_d{{2, 0, 0, 1}};              ///< [I*L^2] current integrate area
constexpr dimension7 magnetic_field_strength_d{{-1, 0, 0, 1}};     ///< [I/L] magnetic moment per volume
constexpr dimension7 magnetization_d = magnetic_field_strength_d;  ///< [I/L] magnetic moment per volume

constexpr dimension7 electric_potential_d{{2, -3, 1, -1}};      ///< [M*L^2/T^3/I] energy per charge
constexpr dimension7 electric_field_strenth_d{{1, -3, 1, -1}};  ///< [M*L/T^3/I] electric potential per length
constexpr dimension7 electric_resistance_d{{2, -3, 1, -2}};     ///< [M*L^2/T^3/I^2] electric potential versus current
constexpr dimension7 electric_conductance_d{{-2, 3, -1, 2}};    ///< [M^-1*L^-2*T^3*I^2] = 1 / electric_resistance
constexpr dimension7 electric_resistivity_d{{3, -3, 1, -2}};    ///< [M*L^3/T^3/I^2] electric_resistance time length
constexpr dimension7 electric_conductivity_d{{-3, 3, -1, 2}};   ///< [M^-1*L^-3*T^3*I^2] 1 / electric_resistivity
constexpr dimension7 electric_capacitance_d{{-2, 4, -1, 2}};  ///< [M^-1*L^-2*T^4*I^2] charge versus electric potential
constexpr dimension7 magnetic_flux_d{{2, -2, 1, -1}};         ///< [M*L^2/T^2/I] energy per current = E/I = B*S
constexpr dimension7 magnetic_flux_density_d{
    {0, -2, 1, -1}};                                ///< [M*T^-2/I] B = electric_field_strenth_d versus velocity
constexpr dimension7 inductance_d{{2, -2, 1, -2}};  ///< [M*L^2/T^2/I^2] magnetic flux versus current, L

constexpr dimension7 electric_chargme_mass_ratio_d{{0, 1, -1, 1}};  ///< [M^-1*T*I]

constexpr dimension7 magnetic_permeability_d{{1, -2, 1, -2}};  ///< [M*L/T^2/I^2], mu
constexpr dimension7 permittivity_d{{-3, 4, -1, 2}};           ///< [M^-1*L^-2*T^4*I^2], epsilon
/// @}

/** \name derived (L, T, M, I, Q) dimension7 */
/// @{
constexpr dimension7 inv_temp_d{{0, 0, 0, 0, -1}};
constexpr dimension7 heat_capacity_d{{2, -2, 1, 0, -1}};  ///< [M*L^2/T^2/Q] energy per temperature
constexpr dimension7 entropy_d = heat_capacity_d;         ///< [M*L^2/T^2/Q] energy per temperature
constexpr dimension7 heat_transfer_coefficient_d{
    {0, -3, 1, 0, -1}};                                            ///< [M/T^3/Q] heat_flux_density versus temperature
constexpr dimension7 specific_heat_capacity_d{{2, -2, 0, 0, -1}};  ///< [L^2/T^2/Q] capacity per mass
constexpr dimension7 thermal_conductivity_d{{1, -3, 1, 0, -1}};    ///< [M*L/T^3/Q]
constexpr dimension7 thermal_insulance_d{{0, 3, -1, 0, 1}};        ///< [M^-1*T^3*Q] = 1 / heat_transfer_coefficient
constexpr dimension7 thermal_resistance_d{{-2, 3, -1, 0, 1}};      ///< [M^-1*L^-2*T^3*Q]
constexpr dimension7 thermal_resistivity_d{{-1, 3, -1, 0, 1}};     ///< [M^-1*L^-1*T^3*Q]
/// @}

/** \name derived (L, T, M, I, Q, N) dimension7 */
/// @{
constexpr dimension7 concentration_d{{-3, 0, 0, 0, 0, 1}};    ///< [N/L^3] amount per volume
constexpr dimension7 molar_energy_d{{2, -2, 1, 0, 0, -1}};    ///< [M*L^2/T^2/N] energy per amount
constexpr dimension7 molar_entropy_d{{2, -2, 1, 0, -1, -1}};  ///< [M*L^2/T^2/Q/N] entropy per amount
/// @}

/** \name derived (L, T, M, I, Q, N, J) dimension7 */
/// @{
constexpr dimension7 luminous_flux_d = luminous_intensity_d;  ///< [J]
constexpr dimension7 illuminance_d{{-2, 0, 0, 0, 0, 0, 1}};   ///< [J/L^2] luminous_intensity per area
constexpr dimension7 luminance_d = illuminance_d;             ///< [J/L^2] luminous_intensity per area
/// @}

const std::map<const dimension7, const std::string> description = {
    // base dimension7
    {dimensionless_d, "dimensionless"},
    {length_d, "length, wavelength"},
    {time_d, "time"},
    {mass_d, "mass"},
    {electric_current_d, "electric_current"},
    {thermodynamic_temperature_d, "thermodynamic_temperature"},
    {amount_of_substance_d, "amount_of_substance"},
    {luminous_intensity_d, "luminous_intensity, luminous_flux"},
    // derived (L, T)
    {wave_number_d, "wave_number"},
    {area_d, "area"},
    {volume_d, "volume"},
    {frequency_d, "frequency"},
    {angular_velocity_d, "angular_velocity"},
    {angular_acceleration_d, "angular_acceleration"},
    {speed_d, "speed"},
    {acceleration_d, "acceleration"},
    {jerk_d, "jerk"},
    {jounce_d, "jounce"},
    {crackle_d, "crackle"},
    {pop_d, "pop"},
    {absement_d, "absement"},
    {area_flow_rate_d, "area_flow_rate, kinematic_viscosity, thermal_diffusivity"},
    {volume_flow_rate_d, "volume_flow_rate"},
    {specific_energy_d, "specific_energy"},
    {absorbed_dose_rate_d, "absorbed_dose_rate"},
    {substance_permeability_d, "substance_permeability"},
    // derived (L, T, M)
    {inertia_d, "inertia"},
    {mass_line_density_d, "mass_line_density"},
    {mass_area_density_d, "mass_area_density"},
    {mass_density_d, "mass_density"},
    {specific_volume_d, "specific_volume"},
    {mass_flow_rate_d, "mass_flow_rate"},
    {mass_flow_acceleration_d, "mass_flow_acceleration, surface_tension"},
    {mass_flow_jerk_d, "mass_flow_jerk, power_area_density"},
    {force_d, "force"},
    {momentum_d, "momentum"},
    {energy_d, "energy, moment_of_force, torque"},
    {angular_momentum_d, "angular_momentum, action"},
    {power_d, "power"},
    {energy_density_d, "energy_density, pressure"},
    {energy_line_density_d, "energy_line_density"},
    {power_density_d, "power_density"},
    {dynamic_viscosity_d, "dynamic_viscosity"},
    // derived (L, T, M, I)
    {current_density_d, "current_density"},
    {electric_charge_d, "electric_charge"},
    {electric_charge_density_d, "electric_charge_density"},
    {electric_area_charge_density_d,
     "electric_area_charge_density, electric_flux_density, electric_displacement_field, electric_polarization_field"},
    {electric_line_charge_density_d, "electric_line_charge_density"},
    {electric_dipole_moment_d, "electric_dipole_moment"},
    {magnetic_moment_d, "magnetic_moment"},
    {magnetic_field_strength_d, "magnetic_field_strength, magnetization"},
    {electric_potential_d, "electric_potential"},
    {electric_field_strenth_d, "electric_field_strenth"},
    {electric_resistance_d, "electric_resistance"},
    {electric_conductance_d, "electric_conductance"},
    {electric_resistivity_d, "electric_resistivity"},
    {electric_conductivity_d, "electric_conductivity"},
    {electric_capacitance_d, "electric_capacitance"},
    {magnetic_flux_d, "magnetic_flux"},
    {magnetic_flux_density_d, "magnetic_flux_density"},
    {inductance_d, "inductance"},
    {electric_chargme_mass_ratio_d, "electric_chargme_mass_ratio"},
    {magnetic_permeability_d, "magnetic_permeability"},
    {permittivity_d, "permittivity"},
    // derived (L, T, M, I, K)
    {heat_capacity_d, "heat_capacity, entropy"},
    {heat_transfer_coefficient_d, "heat_transfer_coefficient"},
    {specific_heat_capacity_d, "specific_heat_capacity"},
    {thermal_conductivity_d, "thermal_conductivity"},
    {thermal_insulance_d, "thermal_insulance"},
    {thermal_resistance_d, "thermal_resistance"},
    {thermal_resistivity_d, "thermal_resistivity"},
    // derived (L, T, M, I, K, N)
    {concentration_d, "concentration"},
    {molar_energy_d, "molar_energy"},
    {molar_entropy_d, "molar_entropy"},
    // derived (L, T, M, I, K, N, J)
    {illuminance_d, "illuminance, luminance"},
};


class uval {  // united value class (value with dimension7; only support */^ operations)
   public:
    using value_type = real_precision;

    uval() : dim(), value(1) {}
    explicit constexpr uval(const value_type& val) : dim{}, value{val} {}
    explicit constexpr uval(const dimension7& dim, const value_type& val = 1L) : dim{dim}, value{val} {}

    dimension7 dim;    ///< dimension7
    value_type value;  ///< magnitude
};
inline constexpr uval operator*(const uval& lhs, const uval& rhs) {
    return uval(lhs.dim * rhs.dim, lhs.value * rhs.value);
}
inline constexpr uval operator*(const real_precision& lhs, const uval& rhs) { return uval(lhs) * rhs; }
inline constexpr uval operator*(const uval& lhs, const real_precision& rhs) { return lhs * uval(rhs); }
inline constexpr uval operator/(const uval& lhs, const uval& rhs) {
    return uval(lhs.dim / rhs.dim, lhs.value / rhs.value);
}
inline constexpr uval operator/(const real_precision& lhs, const uval& rhs) { return uval(lhs) / rhs; }
inline constexpr uval operator/(const uval& lhs, const real_precision& rhs) { return lhs / uval(rhs); }
inline const uval     operator+(const uval& lhs, const uval& rhs) {
    if (lhs.dim != rhs.dim) { throw unit_error("add function must fit the same unit"); }
    return uval(lhs.dim, lhs.value + rhs.value);
}
inline const uval operator-(const uval& lhs, const uval& rhs) {
    if (lhs.dim != rhs.dim) { throw unit_error("minus function must fit the same unit"); }
    return uval(lhs.dim, lhs.value - rhs.value);
}
inline uval power(const uval& lhs, const real_precision& index) {
    return uval(lhs.dim.power(index), std::pow(lhs.value, index));
}
inline bool is_same_dimension7(const uval& lhs, const uval& rhs) { return lhs.dim == rhs.dim; }

inline std::string to_string(const uval& u) {
    std::ostringstream os;
    os << std::setiosflags(std::ios::scientific) << std::setprecision(6) << std::setw(12) << u.value;
    return os.str() + " " + u.dim.to_string();
}

/** \name SI base */
/// @{
constexpr uval _base_1(dimensionless_d);               ///< 1
constexpr uval _base_1m(length_d);                     ///< 1 meter
constexpr uval _base_1s(time_d);                       ///< 1 second
constexpr uval _base_1kg(mass_d);                      ///< 1 kilogram
constexpr uval _base_1A(electric_current_d);           ///< 1 ampere
constexpr uval _base_1K(thermodynamic_temperature_d);  ///< 1 kelvins
constexpr uval _base_1mol(amount_of_substance_d);      ///< 1 mole
constexpr uval _base_1cd(luminous_intensity_d);        ///< 1 candela
// derived
constexpr uval _base_1Hz(frequency_d);             ///< 1 hertz
constexpr uval _base_1N(force_d);                  ///< 1 newton
constexpr uval _base_1Pa(pressure_d);              ///< 1 pascal
constexpr uval _base_1J(energy_d);                 ///< 1 joule
constexpr uval _base_1W(power_d);                  ///< 1 watt
constexpr uval _base_1C(electric_charge_d);        ///< 1 comloub
constexpr uval _base_1V(electric_potential_d);     ///< 1 volt
constexpr uval _base_1F(electric_capacitance_d);   ///< 1 faraday
constexpr uval _base_1S(electric_conductance_d);   ///< 1 siemens
constexpr uval _base_1Om(electric_resistance_d);   ///< 1 ohm
constexpr uval _base_1Wb(magnetic_flux_d);         ///< 1 weber
constexpr uval _base_1T(magnetic_flux_density_d);  ///< 1 tesla
constexpr uval _base_1H(inductance_d);             ///< 1 henry
// others
constexpr uval _nostd_1cal = 4.184L * _base_1J;  ///< 1 cal (non-standard)
/// @}

/** \name universal physical constant */
/// @{
constexpr uval G_gravitional_constant(dimension7{{3, -2, -1}}, 6.6740831E-11L);
constexpr uval c_lightspeed(speed_d, 2.997924580E+8L);
constexpr uval ep0_permittivity(permittivity_d, 8.854187817E-12L);          ///< 1/(4*pi*ke)
constexpr uval mu0_permeability(magnetic_permeability_d, 1.256637061E-6L);  ///< 4*pi*ke/c^2
constexpr uval ke_Comloub(dimensionless_d / permittivity_d, 8.9875517873681764E+9L);
constexpr uval R_gas_constant(molar_entropy_d, 8.314459848L);        ///< k * N
constexpr uval Rydberg_constant(wave_number_d, 10973731.56850865L);  ///< me* e ^ 4 / (8 * ep0 ^ 2 * h ^ 3 * c)
constexpr uval Faraday_constant(electric_charge_d / amount_of_substance_d, 96485.3328959L);  ///< e*N
constexpr uval Stefan_constant(dimension7{{0, -3, 1, -4}}, 5.67036713E-8L);  ///< pi^2 kB^4/(60*hb^3*c^2)
constexpr uval muB_magnetic_moment(magnetic_moment_d, 9.27400999457E-24L);   ///< e*hb / (2*me)
constexpr uval muN_magnetic_moment(magnetic_moment_d, 5.05078369931E-27L);   ///< e*hb / (2*mn)
constexpr uval Bohr_length(length_d, 5.291772106712E-11L);                   ///< hb^2/(ke*me*e^2)
constexpr uval h_Planck(action_d, 6.62607004081E-34L);
constexpr uval hb_Planck(action_d, 1.05457180013E-34L);  ///< h/(2*pi)
constexpr uval me_mass(mass_d, 9.1093835611E-31L);
constexpr uval mp_mass(mass_d, 1.67262189821E-27L);
constexpr uval mn_mass(mass_d, 1.67492749804e-27L);
constexpr uval amu_mass(mass_d, 1.66053886E-27L);
constexpr uval e_charge(electric_charge_d, 1.602176620898E-19L);
constexpr uval k_Boltzman(entropy_d, 1.3806490351E-23L);
constexpr uval N_Avagadro(dimensionless_d / amount_of_substance_d, 6.02214085774E+23L);
/// @}

static const std::map<std::string, real_precision> uval_prefix = {
    {"Y", 1e+24L}, {"Z", 1e+21L}, {"E", 1e+18L}, {"P", 1e+15L}, {"T", 1e+12L}, {"G", 1e+9L},
    {"M", 1e+6L},  {"k", 1e+3L},  {"h", 1e+2L},  {"c", 1e-2L},  {"m", 1e-3L},  {"µ", 1e-6L},
    {"n", 1e-9L},  {"p", 1e-12L}, {"f", 1e-15L}, {"a", 1e-18L}, {"z", 1e-21L}, {"y", 1e-24L},
};
static const std::map<std::string, uval> uval_names = {
    {"1", _base_1},
    {"m", _base_1m},
    {"s", _base_1s},
    {"kg", _base_1kg},
    {"A", _base_1A},
    {"K", _base_1K},
    {"mol", _base_1mol},
    {"cd", _base_1cd},
    {"Hz", _base_1Hz},
    {"N", _base_1N},
    {"Pa", _base_1Pa},
    {"J", _base_1J},
    {"W", _base_1W},
    {"C", _base_1C},
    {"V", _base_1V},
    {"F", _base_1F},
    {"S", _base_1S},
    {"Om", _base_1Om},
    {"Wb", _base_1Wb},
    {"T", _base_1T},
    {"H", _base_1H},
    // non-standard
    {"g", 0.001L * _base_1kg},
    {"Angs", 1e-10L * _base_1m},
    {"Bohr", Bohr_length},
    {"amu", amu_mass},
    {"e", e_charge},
    {"eV", e_charge* _base_1V},
    {"Hart", power(ke_Comloub, 2) * me_mass* power(e_charge, 4) / power(hb_Planck, 2)},
    {"auK", power(ke_Comloub, 2) * me_mass* power(e_charge, 4) / power(hb_Planck, 2) / (k_Boltzman)},
    {"cal", _nostd_1cal},
    {"wn", h_Planck* c_lightspeed / (0.01L * _base_1m)},  ///< alias cm^-1
};

namespace inner {
/** \name minimal constexpr functions */
template <typename T>
CONSTEXPR_DECOR T exp_int(int n) {
    long double x = 2.7182818284590452353602874713526624977572470937000L;
    if (n == 0) return 1;
    if (n < 0) x = 1. / x, n = -n;
    long double y = 1.;
    while (n > 1) {
        if (n % 2 == 0) {
            n /= 2;
        } else {
            y *= x;
            n = (n - 1) / 2;
        }
        x *= x;
    }
    return (T)(x * y);
}
template <typename T>
CONSTEXPR_DECOR long double exp(T num) {
    long double sum   = 0L;
    long double term  = 1L;
    long int    num_i = (long int) (num);
    long double pref  = exp_int<long double>(num_i);
    num -= num_i;
    for (unsigned int count = 1; term >= 1e-100 || count < 100; ++count) {
        sum += term;
        term *= ((long double) num) / count;
    }
    return pref * sum;
}
template <typename T>
CONSTEXPR_DECOR long double log(T num) {
    long double log2 = 0.69314718055994530941723212145817656807550013436026L;
    long double sum  = 0L;
    while (num < 1) num *= 2L, sum -= log2;
    while (num > 2) num /= 2L, sum += log2;
    long double term = ((long double) num - 1) / ((long double) num + 1);
    long double mul  = term * term;
    for (unsigned int tmp_odd = 1; term >= 1e-100 || tmp_odd < 100; tmp_odd += 2) {
        sum += 2 * term / tmp_odd;
        term *= mul;
    }
    return sum;
}
template <typename T>
CONSTEXPR_DECOR long double pow(T a, T b) {
    return inner::exp(b * inner::log(a));
}
/// @}

/** \name minimal linalg utils for solving unit systems at a compile-time cost */
/// @{
template <typename T, std::size_t M, std::size_t N>
struct matrix {
    static_assert(M != 0 && N != 0, "matrix requires positive dimension7");
    using value_type                       = T;
    using size_type                        = std::size_t;
    static constexpr size_type row_size    = M;
    static constexpr size_type column_size = N;
    CONSTEXPR_DECOR T* operator[](std::size_t i) { return _data[i]; }
#ifdef FULLY_COMPILE_TIME_OPT
    CONSTTYPE T const* operator[](std::size_t i) const { return _data[i]; }
#endif  // FULLY_COMPILE_TIME_OPT

    T _data[M][N];
};

// gauss-jordan elimination
template <typename T, std::size_t M, std::size_t N>
CONSTTYPE std::tuple<matrix<T, M, N>, std::size_t, T> gauss_jordan_impl(matrix<T, M, N> m, T tolerance) {
    T           det = 1;
    std::size_t i = 0, j = 0;
    while (i < M && j < N) {
        for (std::size_t ip = i + 1; ip < M; ++ip) {  // Choose largest magnitude as pivot
            if (std::abs(m[ip][j]) > std::abs(m[i][j])) {
                for (std::size_t jp = 0; jp < N; ++jp) {
                    T tmp     = m[i][jp];
                    m[i][jp]  = m[ip][jp];
                    m[ip][jp] = tmp;
                }
                det *= -1;
                break;
            }
        }
        if (!(std::abs(m[i][j]) < tolerance)) {
            // normalization
            auto s = m[i][j];
            for (std::size_t jp = 0; jp < N; ++jp) m[i][jp] /= s;
            det /= s;
            // elimination
            for (std::size_t ip = 0; ip < M; ++ip) {
                if (ip == i) continue;
                if (!(std::abs(m[ip][j]) < tolerance)) {
                    auto s = m[ip][j];
                    for (std::size_t jp = 0; jp < N; ++jp) m[ip][jp] -= s * m[i][jp];
                }
            }
            ++i;  // Select next row (i indicates rank)
        }
        ++j;  // skip && continue to the next column
    }
    det = (i == M) ? det : 0;
    return std::make_tuple(m, i, det);
}

template <typename T, std::size_t M>
CONSTTYPE matrix<T, M, M> inverse(matrix<T, M, M> m) {
    matrix<T, M, M>     inv = {};
    matrix<T, M, 2 * M> mI  = {};
    for (std::size_t i = 0; i < M; ++i) {
        std::size_t j = 0;
        for (std::size_t j1 = 0; j1 < M; ++j1, ++j) mI[i][j] = m[i][j1];
        for (std::size_t j2 = 0; j2 < M; ++j2, ++j) mI[i][j] = (i == j2) ? T{1} : T{0};
    }
    auto t = gauss_jordan_impl(mI, phys::math::eps16);

    if (std::get<1>(t) < M) throw unit_error("not full rank in unit systems");
    auto Iinv = std::get<0>(t);
    for (std::size_t i = 0; i < M; ++i) {
        std::size_t j = M;
        for (std::size_t j2 = 0; j2 < M; ++j2, ++j) inv[i][j2] = Iinv[i][j];
    }
    return inv;
}

template <typename T, std::size_t M, std::size_t N, std::size_t P>
CONSTTYPE matrix<T, M, P> matmul(matrix<T, M, N> a, matrix<T, N, P> b) {
    matrix<T, M, P> c = {0};
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            for (std::size_t k = 0; k < P; ++k) { c[i][k] += a[i][j] * b[j][k]; }
        }
    }
    return c;
}

template <typename T, std::size_t M, std::size_t N>
inline std::ostream& operator<<(std::ostream& os, matrix<T, M, N> m) {
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) os << " " << m[i][j];
        os << "\n";
    }
    return os;
}
/// @}
};  // namespace inner


class unitsys {
   public:
    using value_type = real_precision;

    std::array<uval, dimension7_size> units;
    using MatrixType = inner::matrix<value_type, dimension7_size, dimension7_size>;
    using VectorType = inner::matrix<value_type, dimension7_size, 1>;
    MatrixType invA;

    value_type _c, _h, _hb, _ke, _me, _e, _k, _N, _G;  // intrinsic constant

    explicit CONSTEXPR_DECOR unitsys(const uval& u0, const uval& u1, const uval& u2, const uval& u3, const uval& u4,
                                     const uval& u5, const uval& u6)
        : units{u0, u1, u2, u3, u4, u5, u6},
          invA{inner::inverse(MatrixType{
              u0.dim._data[0], u1.dim._data[0], u2.dim._data[0], u3.dim._data[0],
              u4.dim._data[0], u5.dim._data[0], u6.dim._data[0],  // L
              u0.dim._data[1], u1.dim._data[1], u2.dim._data[1], u3.dim._data[1],
              u4.dim._data[1], u5.dim._data[1], u6.dim._data[1],  // T
              u0.dim._data[2], u1.dim._data[2], u2.dim._data[2], u3.dim._data[2],
              u4.dim._data[2], u5.dim._data[2], u6.dim._data[2],  // M
              u0.dim._data[3], u1.dim._data[3], u2.dim._data[3], u3.dim._data[3],
              u4.dim._data[3], u5.dim._data[3], u6.dim._data[3],  // I
              u0.dim._data[4], u1.dim._data[4], u2.dim._data[4], u3.dim._data[4],
              u4.dim._data[4], u5.dim._data[4], u6.dim._data[4],  // Q
              u0.dim._data[5], u1.dim._data[5], u2.dim._data[5], u3.dim._data[5],
              u4.dim._data[5], u5.dim._data[5], u6.dim._data[5],  // N
              u0.dim._data[6], u1.dim._data[6], u2.dim._data[6], u3.dim._data[6],
              u4.dim._data[6], u5.dim._data[6], u6.dim._data[6]  // J
          })},
          _c{conv(c_lightspeed)},
          _h{conv(h_Planck)},
          _hb{conv(hb_Planck)},
          _ke{conv(ke_Comloub)},
          _me{conv(me_mass)},
          _e{conv(e_charge)},
          _k{conv(k_Boltzman)},
          _N{conv(N_Avagadro)},
          _G{conv(G_gravitional_constant)} {};

    /** \name non-static functions
     * 1) eval a dimension7 in current unitsys
     * 2) conv a uval to current unitsys
     */
    /// @{
    inline CONSTEXPR_DECOR value_type eval(const dimension7 dim) const {
        VectorType vec = {dim._data[0], dim._data[1], dim._data[2], dim._data[3],
                          dim._data[4], dim._data[5], dim._data[6]};
        VectorType sol = inner::matmul(invA, vec);
        value_type val = 1.0L;
        for (int i = 0; i < dimension7_size; ++i) val *= inner::pow(units[i].value, sol._data[i][0]);
        return val;
    }
    inline CONSTEXPR_DECOR value_type conv(const uval& u) { return u.value / this->eval(u.dim); }
    /// @}

    /**
     * \name static functions by calling
     * `phys::us::conv(...)`
     */
    /// @{
    static inline CONSTEXPR_DECOR value_type conv(const dimension7 dim, const unitsys& from_us, const unitsys& to_us) {
        return from_us.eval(dim) / to_us.eval(dim);
    }
    static inline CONSTEXPR_DECOR value_type conv(const uval& u, const unitsys& to_us) {
        return u.value / to_us.eval(u.dim);
    }
    static inline CONSTEXPR_DECOR value_type conv(const unitsys& from_us, const uval& u) {
        return from_us.eval(u.dim) / u.value;
    }

    static uval parse(const std::string& str) {
        size_t pos = str.npos;
        if ((pos = str.find_last_of(" ")) != str.npos) {
            value_type v;
            try {
                v = std::stod(str.substr(0, pos));
            } catch (std::runtime_error& e) { throw unit_error(str + ": parse_error in blankspace"); }
            std::string token2 = str.substr(pos + 1);
            return v * parse(token2);
        }
        if ((pos = str.find_last_of("*")) != str.npos) {
            std::string token1 = str.substr(0, pos);
            std::string token2 = str.substr(pos + 1);
            return parse(token1) * parse(token2);
        }
        if ((pos = str.find_last_of("/")) != str.npos) {
            std::string token1 = str.substr(0, pos);
            std::string token2 = str.substr(pos + 1);
            return parse(token1) / parse(token2);
        }
        if ((pos = str.find_last_of("^")) != str.npos) {
            std::string       token1 = str.substr(0, pos);
            std::string       token2 = str.substr(pos + 1);
            std::stringstream sstr(token2);
            value_type        index;
            sstr >> index;
            return power(parse(token1), index);
        }
        uval u;
        try {
            u = uval_names.at(str);
        } catch (std::out_of_range& e) {
            try {
                value_type v = uval_prefix.at(str.substr(0, 1));
                u            = uval_names.at(str.substr(1));
                u            = u * v;
            } catch (std::out_of_range& e) { throw unit_error(str + ": parse_error in uval"); }
        }
        return u;
    };
    /**
     * Coverting string-style unit [dimension7 used] to a unit system
     */
    static inline value_type conv(const std::string& str, const unitsys& to_us) { return conv(parse(str), to_us); }
    /**
     * Coverting unit system to a string-style unit [dimension7 used]
     */
    static inline value_type conv(const unitsys& from_us, const std::string& str) { return conv(from_us, parse(str)); }
    /**
     * convert different quantities based on energy equivalence priciples
     * [E] = [M]* c^2 = h / [T] = h * c / [L] = kB * [Q] //! not hb but h
     */
    static value_type as(const dimension7 dim, const uval& u, const unitsys& us) {
        value_type ovlp = reduce_l_energy(dim) * reduce_l_energy(u.dim);
        if (ovlp == 0) throw unit_error("invalid conversion");
        uval       ux = (ovlp > 0) ? u : 1L / u;
        value_type vx = unitsys::conv(ux, us);

        value_type diffL = dim._data[0] - ux.dim._data[0];
        value_type diffT = dim._data[1] - ux.dim._data[1];
        value_type diffM = dim._data[2] - ux.dim._data[2];
        value_type diffQ = dim._data[4] - ux.dim._data[4];
        value_type diffN = dim._data[5] - ux.dim._data[5];
        if (dim._data[3] != u.dim._data[3] || dim._data[6] != u.dim._data[6] || diffL + diffT != diffM + diffQ) {
            throw unit_error("invalid conversion");
        }
        value_type coeff_N = -diffN;
        value_type coeff_k = -diffQ;
        value_type coeff_h = diffM - coeff_k;
        value_type coeff_c = -diffT - coeff_h - 2 * coeff_k;

        return vx * std::pow(us._h, coeff_h) * std::pow(us._c, coeff_c) * std::pow(us._k, coeff_k) *
               std::pow(us._N, coeff_N);
    }
    static inline value_type as(const dimension7 dim, const std::string& str, const unitsys& us) {
        return as(dim, parse(str), us);
    }
    /// @}
};

#define PHYS_DEFINE_UNITSYS_NAMESPACE(USNAME, _0, _1, _2, _3, _4, _5, _6)                                         \
                                                                                                                  \
    namespace USNAME {                                                                                            \
    using value_type = real_precision;                                                                            \
    CONSTTYPE unitsys    unit(_0, _1, _2, _3, _4, _5, _6);                                                        \
    CONSTTYPE value_type c  = unitsys::conv(phys::c_lightspeed, unit);                                            \
    CONSTTYPE value_type h  = unitsys::conv(phys::h_Planck, unit);                                                \
    CONSTTYPE value_type hb = unitsys::conv(phys::hb_Planck, unit);                                               \
    CONSTTYPE value_type ke = unitsys::conv(phys::ke_Comloub, unit);                                              \
    CONSTTYPE value_type me = unitsys::conv(phys::me_mass, unit);                                                 \
    CONSTTYPE value_type e  = unitsys::conv(phys::e_charge, unit);                                                \
    CONSTTYPE value_type k  = unitsys::conv(phys::k_Boltzman, unit);                                              \
    CONSTTYPE value_type N  = unitsys::conv(phys::N_Avagadro, unit);                                              \
    CONSTTYPE value_type G  = unitsys::conv(phys::G_gravitional_constant, unit);                                  \
    inline value_type    as(const dimension7 dim, const uval& u) { return unitsys::as(dim, u, unit); }            \
    inline value_type    as(const dimension7 dim, const std::string& str) { return unitsys::as(dim, str, unit); } \
    };

/** \name universal namespace definitions for useful uval systems */
/// @{
PHYS_DEFINE_UNITSYS_NAMESPACE(si, _base_1m, _base_1s, _base_1kg, _base_1A, _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(planck,                 /*planck uval*/
                              c_lightspeed,           /**!< lightspeed */
                              hb_Planck,              /**!< Planck constant */
                              G_gravitional_constant, /**!< Newton's graviton constant */
                              k_Boltzman,             /**!< Boltzman constant */
                              _base_1A, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(god_given,              /*god-given uval*/
                              c_lightspeed,           /**!< lightspeed */
                              hb_Planck,              /**!< Planck constant */
                              G_gravitional_constant, /**!< Newton's graviton constant */
                              _base_1A, _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(stoney,                 /*stoney uval*/
                              c_lightspeed,           /**!< lightspeed */
                              G_gravitional_constant, /**!< Newton's graviton constant */
                              ke_Comloub,             /**!< Comloub constant */
                              e_charge,               /**!< element charge */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(natural,          /*natural uval*/
                              c_lightspeed,     /**!< lightspeed */
                              hb_Planck,        /**!< Planck constant */
                              me_mass,          /**!< mass of electron */
                              ep0_permittivity, /**!< permittivity */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(gauss,              /*gauss uval*/
                              0.01L * _base_1m,   /**!< 1cm */
                              _base_1s,           /**!< 1s */
                              0.001L * _base_1kg, /**!< 1g */
                              ke_Comloub,         /**!< Comloub constant */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(rydberg,                /*rydberg uval*/
                              hb_Planck,              /**!< Planck constant */
                              2 * me_mass,            /**!< mass of electron */
                              e_charge* e_charge / 2, /**!< element charge */
                              ke_Comloub,             /**!< Comloub constant */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(hartree,    /*hartree uval*/
                              hb_Planck,  /**!< Planck constant */
                              me_mass,    /**!< mass of electron */
                              e_charge,   /**!< element charge */
                              ke_Comloub, /**!< Comloub constant */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(ghartree,   /*generalized hartree uval*/
                              hb_Planck,  /**!< Planck constant */
                              me_mass,    /**!< mass of electron */
                              e_charge,   /**!< element charge */
                              ke_Comloub, /**!< Comloub constant */
                              k_Boltzman, /**!< Boltzman constant */
                              N_Avagadro, /**!< Avagadro constant */
                              _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(quantum_chromo_dynamics, /*qcd uval*/
                              c_lightspeed,            /**!< lightspeed */
                              hb_Planck,               /**!< Planck constant */
                              mp_mass,                 /**!< mass of proton */
                              e_charge,                /**!< element charge */
                              _base_1K, _base_1mol, _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(au_test,    /*atomic uval*/
                              hb_Planck,  /**!< Planck constant */
                              me_mass,    /**!< mass of electron */
                              e_charge,   /**!< element charge*/
                              ke_Comloub, /**!< Comloub constant */
                              k_Boltzman, /**!< Boltzman constant */
                              N_Avagadro, /**!< Avagadro constant */
                              _base_1cd);

PHYS_DEFINE_UNITSYS_NAMESPACE(amu,        /*atomic mass uval*/
                              hb_Planck,  /**!< Planck constant */
                              me_mass,    /**!< amu mass */
                              e_charge,   /**!< element charge */
                              ke_Comloub, /**!< Comloub constant */
                              _base_1K, _base_1mol, _base_1cd);

namespace plk = planck;
namespace god = god_given;
namespace sty = stoney;
namespace nat = natural;
namespace cgs = gauss;
namespace ryd = rydberg;
namespace hat = hartree;
namespace au  = ghartree;  ///< generalized hatree is same to hatree when out of statistics (k_Boltzman & N_Avagadro)
namespace qcd = quantum_chromo_dynamics;
/// @}

typedef unitsys us;

/**
 * 1mea means we measure a quantity at 1*N level.
 * in si unit: N = N_Avagadro = 6.022...E+23
 * and in generalized hatree unit (au): N = 1 (or set N_Avagadro as 1)
 *
 * conv() function works well for converting [energy_d/amount_of_substance_d] dimension7
 */
static CONSTTYPE real_precision au_2_amu       = unitsys::conv(au::unit, amu_mass);
static CONSTTYPE real_precision au_2_ang       = unitsys::conv(au::unit, 1e-10L * _base_1m);
static CONSTTYPE real_precision au_2_ev        = unitsys::conv(au::unit, e_charge* _base_1V);
static CONSTTYPE real_precision au_2_J_1mea    = unitsys::conv(au::unit, _base_1J / _base_1mol) * au::N;
static CONSTTYPE real_precision au_2_kcal_1mea = unitsys::conv(au::unit, 1e+3L * _nostd_1cal / _base_1mol) * au::N;
static CONSTTYPE real_precision au_2_g_1mea    = unitsys::conv(au::unit, 1e-3L * _base_1kg / _base_1mol) * au::N;
static CONSTTYPE real_precision au_2_wn        = unitsys::conv(au::unit, h_Planck* c_lightspeed / (0.01L * _base_1m));
static CONSTTYPE real_precision au_2_fs        = unitsys::conv(au::unit, 1e-15L * _base_1s);
static CONSTTYPE real_precision au_2_ps        = unitsys::conv(au::unit, 1e-12L * _base_1s);
static CONSTTYPE real_precision au_2_K         = unitsys::conv(au::unit, _base_1K);
static CONSTTYPE real_precision au_2_angoverps = au_2_ang / au_2_ps;
};  // namespace phys


#undef CONSTTYPE
#undef CONSTEXPR_DECOR
#endif  // PHYS_H

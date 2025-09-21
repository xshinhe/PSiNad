#define CATCH_CONFIG_MAIN
#include <array>
#include <cmath>

#include "catch2/catch_all.hpp"
#include "psnd/phys.h"


using namespace phys;
using namespace phys::math;
using namespace Catch;

// Test suite for mathematical constants in phys::math
TEST_CASE("Mathematical constants in phys::math are correctly defined", "[phys][math][constants]") {
    SECTION("Euler's constant (eu)") { REQUIRE(eu == Approx(2.718281828459045235360287)); }

    SECTION("Pi and its variants") {
        REQUIRE(pi == Approx(3.141592653589793238462643L));
        REQUIRE(twopi == Approx(2 * pi));
        REQUIRE(halfpi == Approx(pi / 2));
    }

    SECTION("Precision epsilon values") {
        REQUIRE(eps8 == Approx(1.0E-8L));
        REQUIRE(eps16 == Approx(1.0E-16L));
        REQUIRE(eps32 == Approx(1.0E-32L));
    }

    SECTION("Square root constants") {
        REQUIRE(sqrttwo == Approx(std::sqrt(2.0)));
        REQUIRE(sqrthalf == Approx(1.0 / std::sqrt(2.0)));
    }

    SECTION("Complex constants") {
        REQUIRE(im == std::complex<real_precision>(0.0L, 1.0L));
        REQUIRE(iu == std::complex<real_precision>(1.0L, 0.0L));
        REQUIRE(iz == std::complex<real_precision>(0.0L, 0.0L));
    }
}

// Test suite for dimensions template class
TEST_CASE("dimensions template class operations", "[phys][dimensions]") {
    using Dim3 = dimensions<int, 3>;
    const Dim3 a{{1, 2, 3}};
    const Dim3 b{{4, 5, 6}};

    SECTION("Equality and inequality operators") {
        REQUIRE(a == a);
        REQUIRE(a != b);
        REQUIRE(Dim3{{1, 2, 3}} == a);
    }

    SECTION("Arithmetic operators") {
        // Test addition
        Dim3 sum = a * b;  // operator* uses array_add
        REQUIRE(sum._data == std::array<int, 3>{{5, 7, 9}});

        // Test subtraction
        Dim3 diff = b / a;  // operator/ uses array_minus
        REQUIRE(diff._data == std::array<int, 3>{{3, 3, 3}});

        // Test scalar multiplication (power)
        Dim3 scaled = a.power(2);  // power uses array_scale
        REQUIRE(scaled._data == std::array<int, 3>{{2, 4, 6}});
    }

    SECTION("Comparison operators") {
        REQUIRE(a < b);
        REQUIRE(b > a);
        REQUIRE(!(a > b));
        REQUIRE(!(b < a));
    }

    SECTION("to_string() method") {
        std::string str = a.to_string();
        REQUIRE(str == "<1, 2, 3, >");  // Matches format: "<i, j, k, >"
    }
}

// Test suite for dimension7 (7-dimensional physical units)
TEST_CASE("dimension7 predefined physical dimensions", "[phys][dimension7]") {
    SECTION("Base dimensions") {
        REQUIRE(dimensionless_d._data == std::array<real_precision, 7>{{0, 0, 0, 0, 0, 0, 0}});
        REQUIRE(length_d._data == std::array<real_precision, 7>{{1, 0, 0, 0, 0, 0, 0}});
        REQUIRE(time_d._data == std::array<real_precision, 7>{{0, 1, 0, 0, 0, 0, 0}});
        REQUIRE(mass_d._data == std::array<real_precision, 7>{{0, 0, 1, 0, 0, 0, 0}});
        REQUIRE(electric_current_d._data == std::array<real_precision, 7>{{0, 0, 0, 1, 0, 0, 0}});
        REQUIRE(thermodynamic_temperature_d._data == std::array<real_precision, 7>{{0, 0, 0, 0, 1, 0, 0}});
        REQUIRE(amount_of_substance_d._data == std::array<real_precision, 7>{{0, 0, 0, 0, 0, 1, 0}});
        REQUIRE(luminous_intensity_d._data == std::array<real_precision, 7>{{0, 0, 0, 0, 0, 0, 1}});
    }

    SECTION("Derived dimensions (L, T)") {
        REQUIRE(wave_number_d._data == std::array<real_precision, 7>{{-1, 0, 0, 0, 0, 0, 0}});
        REQUIRE(area_d._data == std::array<real_precision, 7>{{2, 0, 0, 0, 0, 0, 0}});
        REQUIRE(volume_d._data == std::array<real_precision, 7>{{3, 0, 0, 0, 0, 0, 0}});
        REQUIRE(frequency_d._data == std::array<real_precision, 7>{{0, -1, 0, 0, 0, 0, 0}});
        REQUIRE(speed_d._data == std::array<real_precision, 7>{{1, -1, 0, 0, 0, 0, 0}});
    }

    SECTION("Derived dimensions (L, T, M)") {
        REQUIRE(force_d._data == std::array<real_precision, 7>{{1, -2, 1, 0, 0, 0, 0}});
        REQUIRE(energy_d._data == std::array<real_precision, 7>{{2, -2, 1, 0, 0, 0, 0}});
        REQUIRE(power_d._data == std::array<real_precision, 7>{{2, -3, 1, 0, 0, 0, 0}});
        REQUIRE(pressure_d._data == std::array<real_precision, 7>{{-1, -2, 1, 0, 0, 0, 0}});
    }

    SECTION("dimension7 helper functions") {
        // Test reduce_l_nonzero: counts non-zero dimensions
        REQUIRE(reduce_l_nonzero(length_d) == 1);  // Only L is non-zero
        REQUIRE(reduce_l_nonzero(force_d) == 3);   // L, T, M are non-zero

        // Test reduce_l_energy: custom formula check
        // Formula: -L - T + M + Q
        REQUIRE(reduce_l_energy(energy_d) == Approx(-2 - (-2) + 1 + 0));  // energy_d is [L^2, T^-2, M^1, ...]
    }
}

// Test suite for dimension7 arithmetic operations
TEST_CASE("dimension7 arithmetic operations", "[phys][dimension7][operators]") {
    SECTION("Multiplication (dimension addition)") {
        dimension7 result = length_d * time_d;  // [L] * [T] = [L^1 T^1]
        REQUIRE(result._data == std::array<real_precision, 7>{{1, 1, 0, 0, 0, 0, 0}});
    }

    SECTION("Division (dimension subtraction)") {
        dimension7 result = speed_d / time_d;  // [L T^-1] / [T] = [L T^-2] (acceleration)
        REQUIRE(result == acceleration_d);
    }

    SECTION("Power (scalar multiplication)") {
        dimension7 result = length_d.power(2);  // [L]^2 = [L^2] (area)
        REQUIRE(result == area_d);
    }
}
#define CATCH_CONFIG_MAIN
#include <complex>
#include <string>

#include "catch2/catch_all.hpp"
#include "psnd/Types.h"

// Forward declarations matching Types.h
namespace PROJECT_NS {
class Param {};
class DataSet {};
}  // namespace PROJECT_NS

using namespace Catch;
using namespace PROJECT_NS;

/**
 * @brief Unit tests for type enumeration mapping (as_enum template)
 * Verifies that C++ types are correctly mapped to their psnd_dtype enumeration equivalents
 */
TEST_CASE("as_enum correctly maps types to psnd_dtype", "[Types][as_enum]") {
    CHECK(as_enum<psnd_void>() == psnd_dtype::psnd_void_type);
    CHECK(as_enum<psnd_bool>() == psnd_dtype::psnd_bool_type);
    CHECK(as_enum<psnd_int>() == psnd_dtype::psnd_int_type);
    CHECK(as_enum<psnd_real>() == psnd_dtype::psnd_real_type);
    CHECK(as_enum<psnd_complex>() == psnd_dtype::psnd_complex_type);
    CHECK(as_enum<psnd_str>() == psnd_dtype::psnd_str_type);
    CHECK(as_enum<psnd_param>() == psnd_dtype::psnd_param_type);
    CHECK(as_enum<psnd_dataset>() == psnd_dtype::psnd_dataset_type);
}

/**
 * @brief Unit tests for type string representation (as_str template)
 * Verifies that C++ types are converted to their correct string names
 */
TEST_CASE("as_str returns correct type names", "[Types][as_str]") {
    CHECK(as_str<psnd_void>() == "psnd_void");
    CHECK(as_str<psnd_bool>() == "psnd_bool");
    CHECK(as_str<psnd_int>() == "psnd_int");
    CHECK(as_str<psnd_real>() == "psnd_real");
    CHECK(as_str<psnd_complex>() == "psnd_complex");
    CHECK(as_str<psnd_str>() == "psnd_str");
    CHECK(as_str<psnd_param>() == "psnd_param");
    CHECK(as_str<psnd_dataset>() == "psnd_dataset");
}

/**
 * @brief Unit tests for enum-to-string conversion (enum_t_as_str)
 * Verifies that psnd_dtype enumerators are converted to their correct string names
 */
TEST_CASE("enum_t_as_str returns correct enum names", "[Types][enum_t_as_str]") {
    CHECK(enum_t_as_str(psnd_dtype::psnd_void_type) == "psnd_void");
    CHECK(enum_t_as_str(psnd_dtype::psnd_bool_type) == "psnd_bool");
    CHECK(enum_t_as_str(psnd_dtype::psnd_int_type) == "psnd_int");
    CHECK(enum_t_as_str(psnd_dtype::psnd_real_type) == "psnd_real");
    CHECK(enum_t_as_str(psnd_dtype::psnd_complex_type) == "psnd_complex");
    CHECK(enum_t_as_str(psnd_dtype::psnd_str_type) == "psnd_str");
    CHECK(enum_t_as_str(psnd_dtype::psnd_param_type) == "psnd_param");
    CHECK(enum_t_as_str(psnd_dtype::psnd_dataset_type) == "psnd_dataset");
    CHECK(enum_t_as_str(static_cast<psnd_dtype>(99)) == "unknown");  // Test invalid enum
}

/**
 * @brief Unit tests for type casting utilities
 * Verifies correct behavior of cast, cast_from_complex, and cast_at functions
 */
TEST_CASE("Type casting utilities work correctly", "[Types][casting]") {
    // Test basic numeric casts
    SECTION("Basic type casting") {
        CHECK(cast<int>(5.6) == 5);
        CHECK(cast<psnd_real>(10) == 10.0);
        CHECK(cast<psnd_bool>(0) == false);
        CHECK(cast<psnd_bool>(1) == true);
    }

    // Test complex number casting
    SECTION("Complex number casting") {
        const psnd_complex c(3.0, 4.0);  // (real=3, imag=4)

        // Cast complex to real types (should return real part)
        CHECK(cast<psnd_real>(c) == 3.0);
        CHECK(cast<psnd_int>(c) == 3);

        // Cast complex to complex (identity)
        CHECK(cast<psnd_complex>(c) == c);

        // Explicit cast_from_complex tests
        CHECK(cast_from_complex<psnd_real>(c) == 3.0);
        CHECK(cast_from_complex<psnd_complex>(c) == c);
    }

    // Test pointer-based casting with cast_at
    SECTION("Pointer-based casting (cast_at)") {
        psnd_int int_data[] = {10, 20, 30};
        CHECK(cast_at<psnd_real, psnd_int>(int_data, 1) == 20.0);

        psnd_complex complex_data[] = {psnd_complex(1.5, -2.5), psnd_complex(3.0, 4.0)};
        CHECK(cast_at<psnd_real, psnd_complex>(complex_data, 0) == 1.5);
        CHECK(cast_at<psnd_int, psnd_complex>(complex_data, 1) == 3);
    }
}
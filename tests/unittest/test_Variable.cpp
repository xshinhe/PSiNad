#define CATCH_CONFIG_MAIN
#include <complex>
#include <string>

#include "catch2/catch_all.hpp"
#include "psnd/Shape.h"
#include "psnd/Variable.h"

using namespace Catch;
using namespace PROJECT_NS;

/**
 * @brief Test the string substitution utility function _subreplace
 * Verifies correct replacement of substrings in a source string
 */
TEST_CASE("_subreplace correctly replaces substrings", "[Variable][utility]") {
    SECTION("Replace single occurrence") {
        std::string input  = "a::b::c";
        std::string result = _subreplace(input, "::", ".");
        CHECK(result == "a.b.c");
    }

    SECTION("No occurrence to replace") {
        std::string input  = "abcdef";
        std::string result = _subreplace(input, "x", "y");
        CHECK(result == input);
    }

    SECTION("Empty source string") {
        std::string input  = "";
        std::string result = _subreplace(input, "::", ".");
        CHECK(result.empty());
    }
}

/**
 * @brief Test basic functionality of VARIABLE template class
 * Verifies name, documentation, and shape accessors
 */
TEST_CASE("VARIABLE basic properties", "[Variable][core]") {
    // Create test Shape objects
    Shape static_shape({2, 3});

    // Use Dimension variables for dynamic shape (matches project's actual usage)
    // Avoid manual pointer management - these are managed globally
    std::size_t dim1 = 4, dim2 = 5;
    Shape       dynamic_shape(Shape::dynamic_t{&dim1, &dim2});  // Use dynamic constructor with pointers to local vars

    SECTION("Static shape variable") {
        VARIABLE<int> var("test::static_var", &static_shape, "Test static variable");

        CHECK(var.name() == "test.static_var");      // Verify name substitution
        CHECK(var.doc() == "Test static variable");  // Verify documentation
        CHECK(var.shape().to_string() == "<2,3,>");  // Verify shape access
    }

    SECTION("Dynamic shape variable") {
        VARIABLE<double> var("test::dynamic_var", &dynamic_shape, "Test dynamic variable");

        CHECK(var.name() == "test.dynamic_var");
        CHECK(var.doc() == "Test dynamic variable");
        CHECK(var.shape().to_string() == "<4,5,>");  // Assumes to_string() works as fixed earlier
    }

    // NO manual deletion of pointers - dim1/dim2 are stack variables, not heap-allocated
}

/**
 * @brief Test VARIABLE_BASE::_LIST collection of variables
 * Verifies that all VARIABLE instances are properly registered in the global list
 */
TEST_CASE("VARIABLE_BASE::_LIST collects all variables", "[Variable][registry]") {
    // Store initial list size to restore later
    size_t initial_size = VARIABLE_BASE::_LIST.size();

    // Create temporary variables of different types
    Shape                          s1(1);
    Shape                          s2({2, 2});
    VARIABLE<int>                  var1("list::var1", &s1, "First list variable");
    VARIABLE<std::complex<double>> var2("list::var2", &s2, "Second list variable");

    // Check list growth and content
    CHECK(VARIABLE_BASE::_LIST.size() == initial_size + 2);
    CHECK(VARIABLE_BASE::_LIST.back() == &var2);         // Last added should be var2
    CHECK(VARIABLE_BASE::_LIST[initial_size] == &var1);  // First added should be var1

    // Verify polymorphism through base class pointers
    CHECK(VARIABLE_BASE::_LIST[initial_size]->name() == "list.var1");
    CHECK(VARIABLE_BASE::_LIST[initial_size + 1]->doc() == "Second list variable");

    // Restore initial list state (critical for test isolation)
    while (VARIABLE_BASE::_LIST.size() > initial_size) { VARIABLE_BASE::_LIST.pop_back(); }
}

/**
 * @brief Test VARIABLE template instantiation with various types
 * Verifies template works with common data types used in chemical/condensed matter calculations
 */
TEST_CASE("VARIABLE template supports multiple types", "[Variable][template]") {
    Shape scalar_shape(1);  // Scalar (rank 0) shape

    SECTION("Integer type") {
        VARIABLE<int> var("type::int", &scalar_shape, "Integer variable");
        CHECK(var.name() == "type.int");
    }

    SECTION("Floating point type") {
        VARIABLE<double> var("type::real", &scalar_shape, "Real (double) variable");
        CHECK(var.name() == "type.real");
    }

    SECTION("Complex type") {
        VARIABLE<std::complex<double>> var("type::complex", &scalar_shape, "Complex variable");
        CHECK(var.name() == "type.complex");
    }

    SECTION("Boolean type") {
        VARIABLE<bool> var("type::bool", &scalar_shape, "Boolean variable");
        CHECK(var.name() == "type.bool");
    }
}
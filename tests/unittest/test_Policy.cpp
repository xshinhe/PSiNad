#define CATCH_CONFIG_MAIN
#include <iostream>
#include <sstream>
#include <vector>

#include "catch2/catch_all.hpp"
#include "psnd/Exception.h"
#include "psnd/Policy.h"


// Define test policies INSIDE PROJECT_NS to match project structure
namespace PROJECT_NS {
DEFINE_POLICY(TestPolicy1, ValueA)
DEFINE_POLICY(TestPolicy2, First, Second, Third)
DEFINE_POLICY(TestPolicy3, Alpha, Beta, Gamma, Delta, Epsilon)  // 5 elements
DEFINE_POLICY(EmptyPolicy)                                      // Empty policy (now supported)
};  // namespace PROJECT_NS

using namespace PROJECT_NS;

/**
 * @brief Helper function to capture stdout for testing _help() output
 * Takes a function pointer to the policy's _help() function (avoids using namespaces as template args)
 */
std::string capture_help_output(void (*help_function)()) {
    std::stringstream buffer;
    std::streambuf*   old = std::cout.rdbuf(buffer.rdbuf());
    help_function();  // Call the policy's _help() via function pointer
    std::cout.rdbuf(old);
    return buffer.str();
}

/**
 * @brief Unit tests for Policy macros and generated functionality
 */
TEST_CASE("Policy macro core functionality", "[Policy][core]") {
    SECTION("Single-element policy") {
        CHECK(TestPolicy1::_dict.at("ValueA") == TestPolicy1::ValueA);
        CHECK(TestPolicy1::_dict.size() == 1);
        CHECK(TestPolicy1::_from("ValueA") == TestPolicy1::ValueA);
    }

    SECTION("Multi-element policy") {
        std::vector<std::pair<std::string, TestPolicy2::_type>> expected = {
            {"First", TestPolicy2::First}, {"Second", TestPolicy2::Second}, {"Third", TestPolicy2::Third}};

        CHECK(TestPolicy2::_dict.size() == expected.size());
        for (const auto& [str, val] : expected) { CHECK(TestPolicy2::_dict.at(str) == val); }

        CHECK(TestPolicy2::_from("First") == TestPolicy2::First);
        CHECK(TestPolicy2::_from("Second") == TestPolicy2::Second);
        CHECK(TestPolicy2::_from("Third") == TestPolicy2::Third);
    }

    SECTION("Policy with 5 elements") {
        CHECK(TestPolicy3::_dict.size() == 5);
        CHECK(TestPolicy3::_from("Alpha") == TestPolicy3::Alpha);
        CHECK(TestPolicy3::_from("Epsilon") == TestPolicy3::Epsilon);
    }
}

TEST_CASE("Policy string conversion and error handling", "[Policy][conversion]") {
    SECTION("Invalid string conversion") {
        const std::string invalid_key = "NonExistent";

        // Capture help output using function pointer to TestPolicy2::_help
        std::string help_output = capture_help_output(&TestPolicy2::_help);
        CHECK(help_output.find("Helps for TestPolicy2:") != std::string::npos);
        CHECK(help_output.find("First [available]") != std::string::npos);
        CHECK(help_output.find("Third [available]") != std::string::npos);

        CHECK(TestPolicy2::_from(invalid_key) == TestPolicy2::_type(0));
    }

    SECTION("Case sensitivity (exact match required)") {
        CHECK(TestPolicy2::_from("first") == TestPolicy2::_type(0));   // Lowercase mismatch
        CHECK(TestPolicy2::_from("SECOND") == TestPolicy2::_type(0));  // Uppercase mismatch
    }
}

TEST_CASE("Policy help function output", "[Policy][utility]") {
    SECTION("Help text structure") {
        // Capture help output using function pointer to TestPolicy3::_help
        std::string output = capture_help_output(&TestPolicy3::_help);

        CHECK(output.find("Helps for TestPolicy3:") != std::string::npos);
        CHECK(output.find("Alpha [available]") != std::string::npos);
        CHECK(output.find("Epsilon [available]") != std::string::npos);
    }
}

TEST_CASE("Edge cases for policy macros", "[Policy][edge]") {
    SECTION("Empty policy (no enumerators)") {
        CHECK(EmptyPolicy::_dict.empty() == true);

        // Capture help output using function pointer to EmptyPolicy::_help
        std::string output = capture_help_output(&EmptyPolicy::_help);
        CHECK(output.find("Helps for EmptyPolicy:") != std::string::npos);
        CHECK(EmptyPolicy::_from("anything") == EmptyPolicy::_type(0));
    }
}
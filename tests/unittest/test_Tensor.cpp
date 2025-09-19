#define CATCH_CONFIG_MAIN
#include <complex>
#include <sstream>

#include "catch2/catch_all.hpp"
#include "psnd/Shape.h"
#include "psnd/Tensor.h"
#include "psnd/Types.h"

using namespace Catch;
using namespace PROJECT_NS;

// Helper function to extract numeric values from a string stream
template <typename T>
std::vector<T> extract_values(const std::string& s) {
    std::vector<T>     values;
    std::istringstream iss(s);
    std::string        token;
    bool               found_data = false;

    // Skip type and shape information to find data values
    while (iss >> token) {
        if (token.find_first_of("0123456789.-+") != std::string::npos) {
            found_data = true;
            try {
                if constexpr (std::is_same_v<T, std::complex<double>>) {
                    double real, imag;
                    if (token == "(") {
                        iss >> real;
                        iss >> token;  // Skip comma
                        iss >> imag;
                        iss >> token;  // Skip ")"
                        values.emplace_back(real, imag);
                    }
                } else {
                    values.push_back(static_cast<T>(std::stod(token)));
                }
            } catch (...) {
                // Ignore non-numeric tokens after finding data start
                if (found_data) break;
            }
        }
    }
    return values;
}

TEMPLATE_TEST_CASE("Tensor Basic Initialization", "[Tensor][core][template]", int, double, std::complex<double>) {
    // Test 1D shape
    Shape            shape1d({5});
    Tensor<TestType> tensor1d(shape1d, "1D tensor for testing");

    SECTION("Correctly initializes size from shape") {
        REQUIRE(tensor1d.size() == 5);
        REQUIRE(tensor1d.shape().size() == 5);
        REQUIRE(tensor1d.shape().dims() == std::vector<size_t>{5});
    }

    SECTION("Initializes data to zero") {
        TestType* data = tensor1d.data();
        for (size_t i = 0; i < tensor1d.size(); ++i) { REQUIRE(data[i] == TestType{0}); }
    }

    SECTION("Stores correct documentation info") { REQUIRE(tensor1d.help("any_name") == "1D tensor for testing"); }

    // Test 2D shape
    Shape            shape2d({2, 3});
    Tensor<TestType> tensor2d(shape2d);

    SECTION("Handles multi-dimensional shapes") {
        REQUIRE(tensor2d.size() == 6);  // 2*3
        REQUIRE(tensor2d.shape().dims() == std::vector<size_t>{2, 3});
    }
}

TEMPLATE_TEST_CASE("Tensor Data Manipulation", "[Tensor][core][template]", int, double, std::complex<double>) {
    Shape            shape({3});
    Tensor<TestType> tensor(shape);
    TestType*        data = tensor.data();

    SECTION("Allows modification through data() pointer") {
        data[0] = TestType{1};
        data[1] = TestType{2};
        data[2] = TestType{3};

        REQUIRE(data[0] == TestType{1});
        REQUIRE(data[1] == TestType{2});
        REQUIRE(data[2] == TestType{3});
    }

    SECTION("Maintains data integrity after modification") {
        std::vector<TestType> expected = {TestType{5}, TestType{10}, TestType{15}};
        std::copy(expected.begin(), expected.end(), data);

        for (size_t i = 0; i < expected.size(); ++i) { REQUIRE(data[i] == expected[i]); }
    }
}

TEMPLATE_TEST_CASE("Tensor repr() Method", "[Tensor][core][template]", int, double, std::complex<double>) {
    Shape            shape({2, 2});
    Tensor<TestType> tensor(shape);
    TestType*        data = tensor.data();

    // Populate with known values
    data[0] = TestType{1};
    data[1] = TestType{2};
    data[2] = TestType{3};
    data[3] = TestType{4};

    std::string        representation = tensor.repr();
    std::istringstream repr_stream(representation);
    std::string        first_line;
    std::getline(repr_stream, first_line);

    SECTION("repr() includes correct type information") {
        std::string expected_type = as_str<TestType>();
        REQUIRE(first_line.find(expected_type) != std::string::npos);
    }

    SECTION("repr() includes correct size and shape") {
        REQUIRE(first_line.find("4") != std::string::npos);    // Total elements
        REQUIRE(first_line.find("2 2") != std::string::npos);  // Shape dimensions
    }

    SECTION("repr() contains all data values") {
        std::vector<TestType> extracted = extract_values<TestType>(representation);
        std::vector<TestType> expected  = {TestType{1}, TestType{2}, TestType{3}, TestType{4}};
        REQUIRE(extracted == expected);
    }
}

TEST_CASE("Tensor Type Specific Checks", "[Tensor][core]") {
    SECTION("int Tensor has correct type enum") {
        Tensor<int> tensor(Shape({1}));
        REQUIRE(tensor.type() == as_enum<int>());
    }

    SECTION("double Tensor has correct type enum") {
        Tensor<double> tensor(Shape({1}));
        REQUIRE(tensor.type() == as_enum<double>());
    }

    SECTION("complex<double> Tensor has correct type enum") {
        Tensor<std::complex<double>> tensor(Shape({1}));
        REQUIRE(tensor.type() == as_enum<std::complex<double>>());
    }
}
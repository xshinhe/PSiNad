#define CATCH_CONFIG_MAIN

#include <cstddef>
#include <vector>

#include "catch2/catch_all.hpp"
#include "psnd/Einsum.h"  // Use the correct header path with namespace

using namespace Catch;
using namespace PROJECT_NS;  // Include the project namespace

// Helper function to generate test data (matching the provided example)
std::vector<int> generate_test_data(std::size_t L) {
    std::vector<int> data(L);
    for (std::size_t z = 0; z < L; ++z) {
        data[z] = static_cast<int>(z % 3) - static_cast<int>((z % 5) * (z % 5)) + static_cast<int>(z % 7);
    }
    return data;
}

TEST_CASE("Einsum operations match expected results", "[einsum][tensor]") {
    const std::size_t L = 2 * 2 * 2 * 2 * 3 * 3 * 3 * 3;  // 1296, matching benchmark
    const auto        A = generate_test_data(L);

    SECTION("Identity operation: 'i'") {
        std::vector<int> res(L);
        // Use PROJECT_NS::einsum without explicit template (let deduction work)
        einsum("i", {A.data()}, {{L}}, res.data(), {L});

        // Verify output matches input
        for (std::size_t i = 0; i < L; ++i) { REQUIRE(res[i] == A[i]); }
    }

    SECTION("Sum all elements: 'i->'") {
        std::vector<int> res(1);
        einsum("i->", {A.data()}, {{L}}, res.data(), {1});
        REQUIRE(res[0] == -2589);  // From numpy benchmark
    }

    SECTION("Multi-dimensional contraction: 'ikkkji->j'") {
        const std::vector<std::size_t> input_shape  = {2, 3, 3, 3, 12, 2};
        const std::vector<std::size_t> output_shape = {12};
        std::vector<int>               res(12);

        einsum("ikkkji->j", {A.data()}, {input_shape}, res.data(), output_shape);

        const std::vector<int> expected = {-25, -9, -9, -8, -11, -22, -8, -15, -5, -10, -21, -5};
        for (std::size_t i = 0; i < 12; ++i) { REQUIRE(res[i] == expected[i]); }
    }

    SECTION("Multi-dimensional contraction: 'ikkkji->ik'") {
        const std::vector<std::size_t> input_shape  = {2, 3, 3, 3, 12, 2};
        const std::vector<std::size_t> output_shape = {2, 3};
        std::vector<int>               res(2 * 3);

        einsum("ikkkji->ik", {A.data()}, {input_shape}, res.data(), output_shape);

        const std::vector<int> expected = {-18, -28, -33, -27, -21, -21};
        for (std::size_t i = 0; i < 6; ++i) { REQUIRE(res[i] == expected[i]); }
    }

    SECTION("Vector dot product: 'i,i'") {
        std::vector<int> res(1);
        einsum("i,i", {A.data(), A.data()}, {{L}, {L}}, res.data(), {1});
        REQUIRE(res[0] == 56293);  // From numpy benchmark
    }

    SECTION("Matrix element-wise product sum: 'ik,ik'") {
        const std::vector<std::size_t> input_shape = {16, 81};
        std::vector<int>               res(1);

        einsum("ik,ik", {A.data(), A.data()}, {input_shape, input_shape}, res.data(), {1});
        REQUIRE(res[0] == 56293);
    }

    SECTION("Matrix multiplication sum: 'ik,ki'") {
        const std::vector<std::size_t> shape1 = {16, 81};
        const std::vector<std::size_t> shape2 = {81, 16};
        std::vector<int>               res(1);

        einsum("ik,ki", {A.data(), A.data()}, {shape1, shape2}, res.data(), {1});
        REQUIRE(res[0] == 52891);
    }

    SECTION("Matrix multiplication: 'ik,kj->ij'") {
        const std::vector<std::size_t> shape1       = {4, 324};
        const std::vector<std::size_t> shape2       = {324, 4};
        const std::vector<std::size_t> output_shape = {4, 4};
        std::vector<int>               res(4 * 4);

        einsum("ik,kj->ij", {A.data(), A.data()}, {shape1, shape2}, res.data(), output_shape);

        const std::vector<int> expected = {-2724, 3071,  8230,  7447, -8259, -3074, 3047,  8535,
                                           6290,  -8694, -3109, 3442, 8340,  6031,  -8557, -2810};
        for (std::size_t i = 0; i < 16; ++i) { REQUIRE(res[i] == expected[i]); }
    }

    SECTION("Three-tensor contraction: 'ik,kj,ljjlll->il'") {
        const std::vector<std::size_t> shape1       = {4, 324};
        const std::vector<std::size_t> shape2       = {324, 4};
        const std::vector<std::size_t> shape3       = {3, 4, 4, 3, 3, 3};
        const std::vector<std::size_t> output_shape = {4, 3};
        std::vector<int>               res(4 * 3);

        einsum("ik,kj,ljjlll->il", {A.data(), A.data(), A.data()}, {shape1, shape2, shape3}, res.data(), output_shape);

        const std::vector<int> expected = {83744, 54125,  79687,  57250,  -22579, -1748,
                                           -9172, -21858, -39479, -39026, 55563,  -10344};
        for (std::size_t i = 0; i < 12; ++i) { REQUIRE(res[i] == expected[i]); }
    }
}

TEST_CASE("Einsum handles edge cases correctly", "[einsum][edge]") {
    SECTION("Single element tensor") {
        const std::vector<int> A = {5};
        std::vector<int>       res(1);

        // Test identity
        einsum("a", {A.data()}, {{1}}, res.data(), {1});
        REQUIRE(res[0] == 5);

        // Test sum
        einsum("a->", {A.data()}, {{1}}, res.data(), {1});
        REQUIRE(res[0] == 5);
    }

    SECTION("Empty tensor (zero size)") {
        const std::vector<int> A;
        std::vector<int>       res;

        // Should not crash
        REQUIRE_NOTHROW(einsum("i", {A.data()}, {{0}}, res.data(), {0}));
    }
}
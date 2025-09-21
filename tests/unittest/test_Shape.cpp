#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"
#include "psnd/Shape.h"

using namespace Catch;
using namespace PROJECT_NS;

// Test suite for Shape class basic functionality
TEST_CASE("Shape class basic operations", "[Shape][basic]") {
    SECTION("Constructor from single size_t (rank-1 static shape)") {
        Shape shape(5);

        REQUIRE(shape.rank() == 1);
        REQUIRE(shape.size() == 5);
        REQUIRE(shape() == 5);  // Test operator()
        REQUIRE(shape.dims() == std::vector<std::size_t>{5});
        REQUIRE(shape.to_string() == "<5,>");
    }

    SECTION("Constructor from std::vector (static shape)") {
        std::vector<std::size_t> dims = {2, 3, 4};
        Shape                    shape(dims);

        REQUIRE(shape.rank() == 3);
        REQUIRE(shape.size() == 24);  // 2*3*4
        REQUIRE(shape.dims() == dims);
        REQUIRE(shape.to_string() == "<2,3,4,>");
    }

    SECTION("Static shape with leading dimensions calculation") {
        // Leading dimensions (ldims) are internal, but size depends on them
        Shape shape({2, 2});
        // ldims should be [2, 1], size = 2*2=4
        REQUIRE(shape.size() == 4);
    }

    SECTION("Static shape modification and static_build() impact") {
        Shape shape({2, 3});
        REQUIRE(shape.size() == 6);

        // Modify dims directly (via reference from dims())
        auto& dims = shape.dims();
        dims[0]    = 4;  // Change first dimension to 4

        // Size remains unchanged until static_build() is called (no auto-update for static shape)
        REQUIRE(shape.size() == 6);

        // After static_build(), update() is triggered, recalculating size
        shape.static_build();
        REQUIRE(shape.size() == 12);  // 4*3=12
        REQUIRE(shape.dims() == std::vector<std::size_t>{4, 3});
    }
}

// Test suite for dynamic shape functionality
TEST_CASE("Shape class dynamic shape operations", "[Shape][dynamic]") {
    SECTION("Constructor from pointers (dynamic shape)") {
        std::size_t               a = 2, b = 3;
        std::vector<std::size_t*> dims_ptr = {&a, &b};
        Shape                     shape(dims_ptr);

        // Initial state check
        REQUIRE(shape.rank() == 2);
        REQUIRE(shape.size() == 6);  // 2*3
        REQUIRE(shape.dims() == std::vector<std::size_t>{2, 3});

        // Modify underlying variables and check update
        a = 4;
        REQUIRE(shape.dims() == std::vector<std::size_t>{4, 3});
        REQUIRE(shape.size() == 12);  // 4*3

        b = 5;
        REQUIRE(shape.dims() == std::vector<std::size_t>{4, 5});
        REQUIRE(shape.size() == 20);  // 4*5
    }

    SECTION("Constructor from pointers (dynamic shape with macro)") {
        std::size_t a = 3, b = 3;
        Shape       shape = ShapeDynamic({&a, &b});

        // Initial state check
        REQUIRE(shape.rank() == 2);
        REQUIRE(shape.size() == 9);  // 3*3
        REQUIRE(shape.dims() == std::vector<std::size_t>{3, 3});

        // Modify underlying variables and check update
        a = 4;
        REQUIRE(shape.dims() == std::vector<std::size_t>{4, 3});
        REQUIRE(shape.size() == 12);  // 4*3

        b = 5;
        REQUIRE(shape.dims() == std::vector<std::size_t>{4, 5});
        REQUIRE(shape.size() == 20);  // 4*5
    }

    SECTION("Dynamic shape to static conversion with static_build()") {
        std::size_t x = 1, y = 10;
        Shape       shape(Shape::dynamic_t{&x, &y});

        // Before static_build: dynamic updates work
        x = 2;
        REQUIRE(shape.size() == 20);  // 2*10

        // After static_build: dynamic updates are disabled
        shape.static_build();
        x = 3;  // This change should not affect the shape
        y = 20;
        REQUIRE(shape.dims() == std::vector<std::size_t>{2, 10});
        REQUIRE(shape.size() == 20);
    }

    SECTION("Dynamic shape with underlying dimension validation (assert)") {
        std::size_t invalid_dim = 0;  // Invalid dimension (0 is not allowed)
        Shape       shape(Shape::dynamic_t{&invalid_dim});

        // Accessing dims()/size()/static_build() should trigger assert due to 0 dimension
        REQUIRE_THROWS_AS(shape.size(), std::invalid_argument);
        REQUIRE_THROWS_AS(shape.static_build(), std::invalid_argument);
    }
}

// Test suite for edge cases and special scenarios
TEST_CASE("Shape class edge cases", "[Shape][edge]") {
    SECTION("Rank-1 dynamic shape with value modification") {
        std::size_t s = 100;
        Shape       shape(Shape::dynamic_t{&s});

        REQUIRE(shape.rank() == 1);
        REQUIRE(shape.size() == 100);

        s = 200;
        REQUIRE(shape.size() == 200);
        REQUIRE(shape.dims() == std::vector<std::size_t>{200});
    }


    SECTION("Shape string representation for various ranks") {
        Shape rank1(1);
        REQUIRE(rank1.to_string() == "<1,>");

        Shape rank2({3, 5});
        REQUIRE(rank2.to_string() == "<3,5,>");

        std::size_t* d1 = new std::size_t(2);
        std::size_t* d2 = new std::size_t(4);
        std::size_t* d3 = new std::size_t(6);
        Shape        rank3(Shape::dynamic_t{d1, d2, d3});

        REQUIRE(rank3.to_string() == "<2,4,6,>");
        // Cleanup dynamic pointers to avoid memory leaks
        delete d1;
        delete d2;
        delete d3;
    }

    SECTION("Consistency between dims() and size() after multiple updates") {
        std::size_t a = 2, b = 2, c = 2;
        Shape       shape(Shape::dynamic_t{&a, &b, &c});
        REQUIRE(shape.size() == 8);  // 2*2*2

        a = 3;
        REQUIRE(shape.size() == 12);  // 3*2*2

        c = 3;
        REQUIRE(shape.size() == 18);  // 3*2*3

        b = 4;
        REQUIRE(shape.size() == 36);  // 3*4*3
    }
}
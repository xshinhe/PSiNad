#define CATCH_CONFIG_MAIN
#include <complex>
#include <sstream>

#include "catch2/catch_all.hpp"
#include "psnd/DataSet.h"
#include "psnd/Shape.h"
#include "psnd/Types.h"

using namespace Catch;
using namespace PROJECT_NS;

// Helper function to create a sample DataSet with basic entries
std::unique_ptr<DataSet> create_sample_dataset() {
    auto ds = std::make_unique<DataSet>();
    ds->def_int("int_val", Shape({2, 2}), "2x2 integer matrix");
    ds->def_real("real_val", Shape({3}), "3-element real vector");
    ds->def_complex("complex_val", Shape({1}), "Single complex number");
    ds->def_int("nested.child", 5, "Nested integer");
    return ds;
}

TEST_CASE("DataSet Initialization", "[DataSet][construction]") {
    SECTION("Default constructor creates empty DataSet") {
        DataSet ds;
        REQUIRE(ds._data != nullptr);
        REQUIRE(ds._data->empty() == true);
    }

    SECTION("Initial state has no keys") {
        DataSet ds;
        REQUIRE(ds.haskey("any_key") == false);
    }
}

TEMPLATE_TEST_CASE("DataSet Variable Definition", "[DataSet][def]", psnd_int, psnd_real, psnd_complex) {
    DataSet           ds;
    const std::string key = "test_var";
    Shape             shape({2, 3});
    const std::string info = "Test variable for template";

    SECTION("Define new variable with shape") {
        TestType* ptr = nullptr;
        if constexpr (std::is_same_v<TestType, psnd_int>) {
            ptr = ds.def_int(key, shape, info);
        } else if constexpr (std::is_same_v<TestType, psnd_real>) {
            ptr = ds.def_real(key, shape, info);
        } else if constexpr (std::is_same_v<TestType, psnd_complex>) {
            ptr = ds.def_complex(key, shape, info);
        }

        REQUIRE(ptr != nullptr);
        REQUIRE(ds.haskey(key) == true);

        auto [dtype, data_ptr, shape_ptr] = ds.obtain(key);
        REQUIRE(data_ptr == ptr);
        REQUIRE(*shape_ptr == shape);
        REQUIRE(dtype == as_enum<TestType>());
    }

    SECTION("Define variable with initial array") {
        std::vector<TestType> init_data(shape.size(), TestType{1});
        TestType*             ptr = nullptr;

        if constexpr (std::is_same_v<TestType, psnd_int>) {
            ptr = ds.def_int(key, init_data.data(), shape);
        } else if constexpr (std::is_same_v<TestType, psnd_real>) {
            ptr = ds.def_real(key, init_data.data(), shape);
        } else if constexpr (std::is_same_v<TestType, psnd_complex>) {
            ptr = ds.def_complex(key, init_data.data(), shape);
        }

        REQUIRE(ptr != nullptr);
        for (size_t i = 0; i < shape.size(); ++i) { REQUIRE(ptr[i] == TestType{1}); }
    }

    SECTION("Define variable by referencing existing key") {
        const std::string ref_key = "ref_var";
        if constexpr (std::is_same_v<TestType, psnd_int>) {
            ds.def_int(ref_key, shape);
            ds.def_int(key, ref_key);
        } else if constexpr (std::is_same_v<TestType, psnd_real>) {
            ds.def_real(ref_key, shape);
            ds.def_real(key, ref_key);
        } else if constexpr (std::is_same_v<TestType, psnd_complex>) {
            ds.def_complex(ref_key, shape);
            ds.def_complex(key, ref_key);
        }

        REQUIRE(ds.haskey(key) == true);
        auto [dtype, data_ptr, shape_ptr]             = ds.obtain(key);
        auto [ref_dtype, ref_data_ptr, ref_shape_ptr] = ds.obtain(ref_key);
        REQUIRE(data_ptr == ref_data_ptr);
        REQUIRE(shape_ptr == ref_shape_ptr);
    }

    SECTION("Re-defining existing key throws or replaces correctly") {
        if constexpr (std::is_same_v<TestType, psnd_real>) {
            ds.def_real(key, shape);
            REQUIRE_NOTHROW(ds.def_real_replace(key, shape));
        } else if constexpr (std::is_same_v<TestType, psnd_complex>) {
            ds.def_complex(key, shape);
            REQUIRE_NOTHROW(ds.def_complex_replace(key, shape));
        } else {
            ds.def_int(key, shape);
            REQUIRE_THROWS(ds.def_int(key, shape));  // Regular def throws on redefinition
        }
    }
}

TEST_CASE("DataSet Key Management", "[DataSet][keys]") {
    DataSet           ds;
    const std::string key = "test_key";
    ds.def_int(key, 1);

    SECTION("haskey returns true for existing keys") { REQUIRE(ds.haskey(key) == true); }

    SECTION("haskey returns false for non-existing keys") { REQUIRE(ds.haskey("non_existing") == false); }

    SECTION("undef removes existing key") {
        ds._undef(key);
        REQUIRE(ds.haskey(key) == false);
    }

    SECTION("undef on non-existing key is safe") { REQUIRE_NOTHROW(ds._undef("non_existing")); }

    SECTION("nested keys are properly recognized") {
        ds.def_real("a.b.c", Shape({5}));
        REQUIRE(ds.haskey("a.b.c") == true);
        REQUIRE(ds.haskey("a.b") == true);  // Parent DataSet should exist
        REQUIRE(ds.haskey("a") == true);    // Grandparent DataSet should exist
    }
}

TEST_CASE("DataSet Nested Access", "[DataSet][nested]") {
    auto ds = create_sample_dataset();

    SECTION("at() returns correct nested DataSet") {
        DataSet* nested_ds = ds->at("a");
        REQUIRE(nested_ds != nullptr);
        REQUIRE(nested_ds->haskey("b") == true);
        REQUIRE(nested_ds->haskey("c") == true);
    }

    SECTION("at() on non-existing nested path returns nullptr") { REQUIRE(ds->at("non.existing.path") == nullptr); }

    SECTION("node() returns valid Node for existing keys") {
        Node* node = ds->node("int_val");
        REQUIRE(node != nullptr);
    }

    SECTION("node() returns nullptr for non-existing keys") { REQUIRE(ds->node("non_existing") == nullptr); }
}

TEST_CASE("DataSet Information Retrieval", "[DataSet][info]") {
    DataSet           ds;
    const std::string key  = "info_test";
    const std::string info = "This is a test variable";
    ds.def_real(key, Shape({1}), info);

    SECTION("help() returns correct information") { REQUIRE(ds.help(key) == info); }

    SECTION("help() on non-existing key returns empty string") { REQUIRE(ds.help("non_existing").empty() == true); }

    SECTION("obtain() returns correct type, data and shape") {
        auto [dtype, data_ptr, shape_ptr] = ds.obtain(key);
        REQUIRE(dtype == psnd_dtype::psnd_real_type);
        REQUIRE(data_ptr != nullptr);
        REQUIRE(*shape_ptr == Shape({1}));
    }

    SECTION("obtain() on non-existing key throws") { REQUIRE_THROWS(ds.obtain("non_existing")); }
}

TEST_CASE("DataSet Serialization", "[DataSet][io]") {
    auto ds = create_sample_dataset();

    SECTION("repr() contains key information") {
        std::string representation = ds->repr();
        REQUIRE(representation.find("int_val") != std::string::npos);
        REQUIRE(representation.find("real_val") != std::string::npos);
        REQUIRE(representation.find("complex_val") != std::string::npos);
        REQUIRE(representation.find("nested.child") != std::string::npos);
    }

    SECTION("dump() writes valid data") {
        std::stringstream ss;
        ds->dump(ss);
        std::string dump_str = ss.str();
        REQUIRE(dump_str.empty() == false);
        REQUIRE(dump_str.find("int_val") != std::string::npos);
    }

    SECTION("load() reconstructs DataSet correctly") {
        std::stringstream ss;
        ds->dump(ss);

        DataSet loaded_ds;
        loaded_ds.load(ss);

        REQUIRE(loaded_ds.haskey("int_val") == true);
        REQUIRE(loaded_ds.haskey("real_val") == true);
        REQUIRE(loaded_ds.haskey("complex_val") == true);
        REQUIRE(loaded_ds.haskey("nested.child") == true);
    }
}

TEST_CASE("DataSet Advanced Operations", "[DataSet][advanced]") {
    DataSet ds;

    SECTION("def() with VARIABLE and span works") {
        Shape                  shape({3});
        VARIABLE<psnd_real>    var("var.span", &shape, "test for span works");
        std::vector<psnd_real> data = {1.0, 2.0, 3.0};
        span<psnd_real>        data_span(data.data(), data.size());

        auto result_span = ds.def(var, data_span);
        REQUIRE(result_span.size() == 3);
        REQUIRE(ds.haskey("var.span") == true);

        auto [dtype, data_ptr, shape_ptr] = ds.obtain("var.span");
        REQUIRE(*static_cast<psnd_real*>(data_ptr) == 1.0);
    }

    SECTION("def() with shared DataSet works") {
        auto sub_ds = std::make_shared<DataSet>();
        sub_ds->def_int("sub.key", 42);

        ds.def(sub_ds);
        REQUIRE(ds.haskey("sub.key") == true);
        REQUIRE(ds.at("sub") != nullptr);
    }

    SECTION("load_reframe() handles resizing") {
        std::stringstream ss;
        auto              original = std::make_unique<DataSet>();
        original->def_real("frame_data", Shape({2}), "Frame data");
        original->dump(ss);

        DataSet reframed;
        reframed.load_reframe(ss, 3);  // Resize to 3 samples
        auto [dtype, data_ptr, shape_ptr] = reframed.obtain("frame_data");
        REQUIRE(shape_ptr->size() == 2 * 3);  // Original size * nsamp
    }
}
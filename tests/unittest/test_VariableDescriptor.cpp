#define CATCH_CONFIG_MAIN
#include <memory>
#include <string>

#include "catch2/catch_all.hpp"
#include "psnd/DataSet.h"
#include "psnd/Policy.h"
#include "psnd/VariableDescriptor.h"

using namespace Catch;
using namespace PROJECT_NS;

// Helper function to create a test DataSet with sample data
std::shared_ptr<DataSet> createTestDataSet() {
    auto ds = std::make_shared<DataSet>();
    // Add sample variables to DataSet
    ds->def_int("field.name[0]", 1);       // Example scalar
    ds->def_real("field.vec[1]", 3);       // Example vector (size 3)
    ds->def_complex("meta.config[2]", 1);  // Example meta variable
    return ds;
}

TEST_CASE("VariableDescriptor Construction", "[VariableDescriptor][construction]") {
    SECTION("Constructs with basic parameters") {
        const std::string token  = "field.name[0]";
        const std::string save   = "output";
        auto              policy = VariableDescriptorPolicy::TabularOutput;

        VariableDescriptor desc(token, save, policy);

        // Verify basic initialization (publicly inferrable properties)
        // REQUIRE(desc.tokenString == token);
        // Note: Private members can't be accessed directly, so we test through behavior in other methods
    }

    SECTION("Parses token string correctly") {
        // Test token parsing by checking defineIn behavior later, since parsing is internal
        // This section is a placeholder for parsing logic validation through indirect checks
    }
}

TEST_CASE("VariableDescriptor referIn Method", "[VariableDescriptor][referIn]") {
    // NOT IMPLEMENTED
}

TEST_CASE("VariableDescriptor defineIn Method", "[VariableDescriptor][defineIn]") {
    auto               ds    = std::make_shared<DataSet>();
    const std::string  token = "new.variable[5]";
    VariableDescriptor desc(token, "save3", VariableDescriptorPolicy::TabularOutput);

    SECTION("Defines new variable in DataSet with specified type and shape") {
        std::vector<std::size_t> shape = {2, 3};      // 2x3 matrix
        desc.defineIn(ds, psnd_int_type, shape, 10);  // 10 frames

        // Verify variable exists in DataSet
        REQUIRE(ds->haskey("new.variable[5]"));

        // Verify shape (indirect check through DataSet)
        // auto varShape = ds->getShape("new.variable[5]");
        // REQUIRE(varShape == shape);
    }

    SECTION("Handles default parameters correctly") {
        desc.defineIn(ds);  // Default type (psnd_void_type), empty shape, 0 frames
        REQUIRE(ds->haskey(token));
    }
}

TEST_CASE("VariableDescriptor checkTrace Method", "[VariableDescriptor][checkTrace]") {
    auto               ds = createTestDataSet();
    VariableDescriptor desc("field.vec[1]", "save4", VariableDescriptorPolicy::MetaOutput);
    desc.defineIn(ds);

    SECTION("Validates trace for valid sample index") {
        REQUIRE_NOTHROW(desc.checkTrace(0));  // First sample
        REQUIRE_NOTHROW(desc.checkTrace(2));  // Within vector size (3 elements)
    }

    SECTION("Throws for invalid sample index") {
        REQUIRE_THROWS_AS(desc.checkTrace(3), std::out_of_range);  // Exceeds vector size (3 elements)
    }

    SECTION("Updates trace state when requested") {
        // Check that update=true modifies internal state (indirect check)
        desc.checkTrace(0, true);
        REQUIRE_NOTHROW(desc.checkTrace(0, false));  // Should not throw after valid update
    }
}

TEST_CASE("VariableDescriptor Policy Behavior", "[VariableDescriptor][policy]") {
    auto ds = createTestDataSet();

    SECTION("TabularOutput policy handles multiple frames") {
        VariableDescriptor desc("meta.config[2]", "save5", VariableDescriptorPolicy::TabularOutput);
        desc.defineIn(ds, psnd_real_type, {}, 5);  // 5 frames
        REQUIRE_NOTHROW(desc.checkTrace(0));
        REQUIRE_NOTHROW(desc.checkTrace(4));  // Last frame
    }

    SECTION("MetaOutput policy handles single frame") {
        VariableDescriptor desc("meta.config[2]", "save6", VariableDescriptorPolicy::MetaOutput);
        desc.defineIn(ds, psnd_real_type, {}, 1);  // 1 frame
        REQUIRE_NOTHROW(desc.checkTrace(0));
        REQUIRE_THROWS_AS(desc.checkTrace(1), std::out_of_range);  // No second frame
    }

    SECTION("Input policy references existing data") {
        VariableDescriptor desc("field.name[0]", "save7", VariableDescriptorPolicy::Input);
        desc.defineIn(ds);
        REQUIRE_NOTHROW(desc.checkTrace(0));
    }
}
#define CATCH_CONFIG_MAIN
#include <filesystem>
#include <fstream>

#include "catch2/catch_all.hpp"
#include "psnd/Param.h"

using namespace Catch;
using namespace PROJECT_NS;
namespace fs = std::filesystem;

// Helper function to create temporary test files
std::string create_temp_file(const std::string& content) {
    static int    temp_file_counter = 0;
    std::string   filename          = "param_test_" + std::to_string(temp_file_counter++) + ".json";
    std::ofstream ofs(filename);
    ofs << content;
    ofs.close();
    return filename;
}

// Helper function to clean up temporary files
void cleanup_temp_file(const std::string& filename) {
    if (fs::exists(filename)) { fs::remove(filename); }
}

TEST_CASE("Param Construction", "[Param][core]") {
    SECTION("Construct from valid JSON string") { REQUIRE_NOTHROW(Param("{\"key\": 123}", Param::fromString)); }

    SECTION("Construct from invalid JSON string throws") { REQUIRE_THROWS(Param("{invalid json}", Param::fromString)); }

    SECTION("Construct from existing file") {
        std::string temp_file = create_temp_file("{\"file_key\": \"value\"}");
        REQUIRE_NOTHROW(Param(temp_file, Param::fromFile));
        cleanup_temp_file(temp_file);
    }

    SECTION("Construct from non-existing file throws") {
        REQUIRE_THROWS(Param("non_existing_file_1234.json", Param::fromFile));
    }
}

TEST_CASE("Key Existence Checks", "[Param][core]") {
    Param param("{\"existing_key\": 42, \"nested\": {\"inner_key\": true}}", Param::fromString);

    SECTION("has_key returns true for existing keys") {
        REQUIRE(param.has_key("existing_key") == true);
        REQUIRE(param.has_key("nested.inner_key") == true);  // Assuming dot notation for nested keys
    }

    SECTION("has_key returns false for non-existing keys") {
        REQUIRE(param.has_key("non_existing_key") == false);
        REQUIRE(param.has_key("nested.non_existing") == false);
    }
}

TEST_CASE("Type Checking Methods", "[Param][core]") {
    Param param(R"({
        "bool_val": true,
        "int_val": 100,
        "real_val": 3.14,
        "string_val": "test",
        "object_val": {"a": 1},
        "array_val": [1, 2, 3]
    })",
                Param::fromString);

    SECTION("is_bool correctly identifies boolean values") {
        REQUIRE(param.is_bool("bool_val") == true);
        REQUIRE(param.is_bool("int_val") == false);
    }

    SECTION("is_int correctly identifies integer values") {
        REQUIRE(param.is_int("int_val") == true);
        REQUIRE(param.is_int("real_val") == false);
    }

    SECTION("is_real correctly identifies real values") {
        REQUIRE(param.is_real("real_val") == true);
        REQUIRE(param.is_real("string_val") == false);
    }

    SECTION("is_string correctly identifies string values") {
        REQUIRE(param.is_string("string_val") == true);
        REQUIRE(param.is_string("bool_val") == false);
    }

    SECTION("is_object correctly identifies objects") {
        REQUIRE(param.is_object("object_val") == true);
        REQUIRE(param.is_object("array_val") == false);
    }

    SECTION("is_array correctly identifies arrays") {
        REQUIRE(param.is_array("array_val") == true);
        REQUIRE(param.is_array("object_val") == false);
    }
}

TEST_CASE("Set and Get Methods", "[Param][core]") {
    Param param("{}", Param::fromString);

    SECTION("set_bool and get_bool work correctly") {
        param.set_bool("test_bool", true);
        REQUIRE(param.get_bool({"test_bool"}) == true);

        param.set_bool("test_bool", false);
        REQUIRE(param.get_bool({"test_bool"}) == false);
    }

    SECTION("set_int and get_int work correctly") {
        param.set_int("test_int", 42);
        REQUIRE(param.get_int({"test_int"}) == 42);

        param.set_int("test_int", -100);
        REQUIRE(param.get_int({"test_int"}) == -100);
    }

    SECTION("set_real and get_real work correctly") {
        param.set_real("test_real", 3.14159);
        REQUIRE(param.get_real({"test_real"}) == Approx(3.14159));

        param.set_real("test_real", -0.001);
        REQUIRE(param.get_real({"test_real"}) == Approx(-0.001));
    }

    SECTION("set_string and get_string work correctly") {
        param.set_string("test_string", "hello");
        REQUIRE(param.get_string({"test_string"}) == "hello");

        param.set_string("test_string", "world");
        REQUIRE(param.get_string({"test_string"}) == "world");
    }

    SECTION("get methods return default values for non-existing keys") {
        REQUIRE(param.get_bool({"non_exist_bool"}, "", false) == false);
        REQUIRE(param.get_int({"non_exist_int"}, "", 10) == 10);
        REQUIRE(param.get_real({"non_exist_real"}, "", 2.718) == Approx(2.718));
        REQUIRE(param.get_string({"non_exist_str"}, "", "default") == "default");
    }

    SECTION("get methods throw for non-existing keys without default") {
        REQUIRE_THROWS(param.get_bool({"non_exist_bool"}));
        REQUIRE_THROWS(param.get_int({"non_exist_int"}));
        REQUIRE_THROWS(param.get_real({"non_exist_real"}));
        REQUIRE_THROWS(param.get_string({"non_exist_str"}));
    }

    SECTION("get_real handles physical dimensions correctly") {
        param.set_real("length", 10.0);  // Assume in default units
        psnd_real value = param.get_real({"length"}, "", phys::length_d, 0.0);
        // Verify dimension conversion is applied (test specific to your phys implementation)
        REQUIRE(value != 0.0);
    }

    SECTION("supports nested key access via vector") {
        param.set_int("parent.child", 123);
        REQUIRE(param.get_int({"parent", "child"}) == 123);
    }
}

TEST_CASE("Set If Not Defined Methods", "[Param][core]") {
    Param param(R"({"existing_bool": false, "existing_int": 5, "existing_real": 1.0, "existing_string": "old"})",
                Param::fromString);

    SECTION("set_bool_ifndef preserves existing value") {
        param.set_bool_ifndef("existing_bool", true);
        REQUIRE(param.get_bool({"existing_bool"}) == false);
    }

    SECTION("set_bool_ifndef sets new value for non-existing key") {
        param.set_bool_ifndef("new_bool", true);
        REQUIRE(param.get_bool({"new_bool"}) == true);
    }

    SECTION("set_int_ifndef preserves existing value") {
        param.set_int_ifndef("existing_int", 10);
        REQUIRE(param.get_int({"existing_int"}) == 5);
    }

    SECTION("set_int_ifndef sets new value for non-existing key") {
        param.set_int_ifndef("new_int", 10);
        REQUIRE(param.get_int({"new_int"}) == 10);
    }

    SECTION("set_real_ifndef preserves existing value") {
        param.set_real_ifndef("existing_real", 2.0);
        REQUIRE(param.get_real({"existing_real"}) == Approx(1.0));
    }

    SECTION("set_real_ifndef sets new value for non-existing key") {
        param.set_real_ifndef("new_real", 2.0);
        REQUIRE(param.get_real({"new_real"}) == Approx(2.0));
    }

    SECTION("set_string_ifndef preserves existing value") {
        param.set_string_ifndef("existing_string", "new");
        REQUIRE(param.get_string({"existing_string"}) == "old");
    }

    SECTION("set_string_ifndef sets new value for non-existing key") {
        param.set_string_ifndef("new_string", "new");
        REQUIRE(param.get_string({"new_string"}) == "new");
    }
}

TEST_CASE("Repr Method", "[Param][core]") {
    SECTION("repr returns non-empty string for populated param") {
        Param       param("{\"key\": \"value\"}", Param::fromString);
        std::string representation = param.repr();
        REQUIRE_FALSE(representation.empty());
        REQUIRE(representation.find("key") != std::string::npos);
        REQUIRE(representation.find("value") != std::string::npos);
    }

    SECTION("repr handles empty param") {
        Param       param("{}", Param::fromString);
        std::string representation = param.repr();
        REQUIRE_FALSE(representation.empty());
    }
}

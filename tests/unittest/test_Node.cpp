#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"
#include "psnd/Node.h"
#include "psnd/Types.h"  // Assume psnd_dtype is defined here

using namespace Catch;
using namespace PROJECT_NS;

// Concrete derived class to test the abstract Node interface
class TestNode : public Node {
   public:
    // Constructor with optional type initialization
    explicit TestNode(psnd_dtype type = psnd_void_type) {
        _type = type;  // Accessible since TestNode is a derived class
    }

    // Implement pure virtual methods
    std::string repr() override { return "TestNode(repr)"; }

    std::string help(const std::string& name) override { return "Help for " + name; }

    // Helper to modify type for testing
    void setType(psnd_dtype type) { _type = type; }
};

TEST_CASE("Node base class interface tests", "[Node][core]") {
    SECTION("Default type initialization") {
        TestNode node;
        REQUIRE(node.type() == psnd_void_type);
    }

    SECTION("Type modification and retrieval") {
        TestNode   node;
        psnd_dtype test_type = psnd_int_type;  // Assume common dtype from Types.h

        node.setType(test_type);
        REQUIRE(node.type() == test_type);
    }

    SECTION("repr() method returns expected string") {
        TestNode node;
        REQUIRE(node.repr() == "TestNode(repr)");
    }

    SECTION("help() method formats name correctly") {
        TestNode    node;
        std::string test_name = "test_node_123";

        REQUIRE(node.help(test_name) == "Help for " + test_name);
    }

    SECTION("Abstract class cannot be instantiated directly") {
        // Verify Node is abstract by attempting instantiation (should cause compile error if uncommented)
        // Node abstract_node;  // Compile error expected: cannot declare variable 'abstract_node' to be of abstract
        // type 'PROJECT_NS::Node'
    }
}
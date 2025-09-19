#define CATCH_CONFIG_MAIN
#include <memory>
#include <string>

#include "catch2/catch_all.hpp"
#include "psnd/DataSet.h"
#include "psnd/Kernel.h"
#include "psnd/Param.h"
#include "psnd/RuleSet.h"
#include "psnd/Status.h"

using namespace Catch;
using namespace PROJECT_NS;

// Helper derived class to access protected members for testing
class TestKernel : public Kernel {
   public:
    using Kernel::_child_kernels;
    using Kernel::_dataset;
    using Kernel::_param;
    using Kernel::_ruleset;
    using Kernel::count_calc;
    using Kernel::count_exec;
    using Kernel::depth;
    using Kernel::enable_call_child;
    using Kernel::exec_time;
    using Kernel::has_parent;
    using Kernel::is_timing;
    using Kernel::Kernel;  // Inherit constructors
    using Kernel::kernel_id;
    using Kernel::kernel_name;
    using Kernel::kernel_type;
    using Kernel::max_align_size;
    using Kernel::once_called;

    // Expose protected virtual methods
    void    setInputParam_impl(std::shared_ptr<Param> PM) override { Kernel::setInputParam_impl(PM); }
    void    setInputDataSet_impl(std::shared_ptr<DataSet> DS) override { Kernel::setInputDataSet_impl(DS); }
    Status& initializeKernel_impl(Status& stat) override { return Kernel::initializeKernel_impl(stat); }
    Status& executeKernel_impl(Status& stat) override { return Kernel::executeKernel_impl(stat); }
    Status& finalizeKernel_impl(Status& stat) override { return Kernel::finalizeKernel_impl(stat); }
};


// Helper to create a test DataSet
std::shared_ptr<DataSet> createTestDataSet() { return std::make_shared<DataSet>(); }

// Helper to create a RuleSet by registering a RuleEvaluator
std::shared_ptr<RuleSet> createTestRuleSet() {
    auto ds = createTestDataSet();
    // Create a dummy RuleEvaluator
    auto exprRule = std::make_shared<RuleEvaluator>("result = var1 + var2",  // Dummy rule expression
                                                    ds,                      // Associated DataSet
                                                    "average",               // Mode
                                                    "test_ruleset.dat",      // Save path
                                                    10                       // Total frames
    );
    // Register the RuleEvaluator to create a RuleSet
    RuleSet::registerRulesInRuleSet(exprRule);
    // Retrieve the newly created RuleSet from the global list
    auto& ruleSets = RuleSet::getRuleSets();
    REQUIRE(!ruleSets.empty());
    return ruleSets.back();  // Return the latest RuleSet
}


TEST_CASE("Kernel RuleSet Association", "[Kernel][RuleSet]") {
    Kernel ker("rule_kernel");
    auto   testRuleSet = createTestRuleSet();

    // Test setting a RuleSet on the kernel
    ker.setRuleSet(testRuleSet);
    REQUIRE(ker.getRuleSet() == testRuleSet);  // Verify association

    // Test RuleSet propagation to child kernels
    auto childKer = std::make_shared<Kernel>("child_kernel");
    ker.appendChild(childKer);
    REQUIRE(childKer->getRuleSet() == testRuleSet);  // Child inherits RuleSet
}


TEST_CASE("Kernel Construction and Basic Properties", "[Kernel][construction]") {
    SECTION("Default constructor initializes with empty name") {
        TestKernel ker;
        REQUIRE(ker.getName() == "Kernel__");
        REQUIRE(ker.kernel_name.empty());
        REQUIRE(ker.kernel_id >= 0);
        REQUIRE(ker.depth == 0);
        REQUIRE(ker._child_kernels.empty());
    }

    SECTION("Custom name constructor sets name correctly") {
        const std::string test_name = "SimulationKernel";
        TestKernel        ker(test_name);
        REQUIRE(ker.getName() == "Kernel__" + test_name);
        REQUIRE(ker.kernel_name == test_name);
    }

    SECTION("Each kernel gets unique ID") {
        TestKernel ker1;
        TestKernel ker2;
        REQUIRE(ker1.kernel_id != ker2.kernel_id);
    }
}

TEST_CASE("Kernel Parameter and DataSet Management", "[Kernel][data]") {
    TestKernel ker("TestKernel");
    auto       param   = std::make_shared<Param>("{}", Param::fromString);
    auto       dataset = std::make_shared<DataSet>();

    SECTION("setInputParam correctly initializes parameters") {
        ker.setInputParam(param);
        REQUIRE(ker.getParam() == param);
        REQUIRE(ker._param == param);
    }

    SECTION("setInputDataSet correctly initializes dataset") {
        ker.setInputDataSet(dataset);
        REQUIRE(ker.getDataSet() == dataset);
        REQUIRE(ker._dataset == dataset);
    }

    SECTION("Initial state has no parameters or dataset") {
        REQUIRE(ker.getParam() == nullptr);
        REQUIRE(ker.getDataSet() == nullptr);
        REQUIRE(ker._param == nullptr);
        REQUIRE(ker._dataset == nullptr);
    }
}

TEST_CASE("Kernel Child Management", "[Kernel][children]") {
    TestKernel parent("Parent");
    auto       child1 = std::make_shared<TestKernel>("Child1");
    auto       child2 = std::make_shared<TestKernel>("Child2");
    auto       child3 = std::make_shared<TestKernel>("Child3");

    SECTION("appendChild adds child to hierarchy") {
        parent.appendChild(child1);
        REQUIRE(parent._child_kernels.size() == 1);
        REQUIRE(parent._child_kernels[0] == child1);
        REQUIRE(child1->has_parent == true);

        parent.appendChild(child2);
        REQUIRE(parent._child_kernels.size() == 2);
        REQUIRE(parent._child_kernels[1] == child2);

        // Verify child count (indirectly via child IDs or internal state)
        // Note: Kernel doesn't expose child count publicly, so we test behavior
        REQUIRE(child1->getID() != child2->getID());  // Unique IDs
        REQUIRE(child1->getID() != parent.getID());   // Parent and child have different IDs
    }

    SECTION("insertAt places child at specified index") {
        parent.appendChild(child1);    // Index 0
        parent.insertAt({0}, child2);  // Insert at 0, child1 moves to 1

        REQUIRE(parent._child_kernels.size() == 2);
        REQUIRE(parent._child_kernels[0] == child2);
        REQUIRE(parent._child_kernels[1] == child1);
        REQUIRE(child2->has_parent == true);
    }

    SECTION("removeAt removes specified child") {
        parent.appendChild(child1);
        parent.appendChild(child2);
        parent.appendChild(child3);

        parent.removeAt({1});  // Remove child2
        REQUIRE(parent._child_kernels.size() == 2);
        REQUIRE(parent._child_kernels[0] == child1);
        REQUIRE(parent._child_kernels[1] == child3);
    }

    SECTION("updateAt replaces child at specified index") {
        parent.appendChild(child1);
        parent.updateAt({0}, child2);

        REQUIRE(parent._child_kernels.size() == 1);
        REQUIRE(parent._child_kernels[0] == child2);
        REQUIRE(child2->has_parent == true);
    }

    SECTION("getLastParentKernelAndChildOrder returns correct parent info") {
        parent.appendChild(child1);
        auto [parent_ptr, order] = child1->getLastParentKernelAndChildOrder();

        REQUIRE(parent_ptr == &parent);
        REQUIRE(order == 0);
    }
}

TEST_CASE("Kernel Lifecycle Methods", "[Kernel][lifecycle]") {
    TestKernel ker("LifecycleKernel");
    Status     stat;
    auto       ds    = createTestDataSet();
    auto       param = std::make_shared<Param>("{}", Param::fromString);

    // Test setting input parameters and data
    ker.setInputParam(param);
    ker.setInputDataSet(ds);
    REQUIRE(ker.getParam() == param);
    REQUIRE(ker.getDataSet() == ds);

    SECTION("initializeKernel increments calculation counter") {
        ker.initializeKernel(stat);
        REQUIRE(ker.count_calc == 1);
        REQUIRE(stat.icalc == 0);  // Default initial value
    }

    SECTION("executeKernel increments execution counter") {
        ker.executeKernel(stat);
        REQUIRE(ker.count_exec == 1);
    }

    SECTION("finalizeKernel completes lifecycle without errors") { REQUIRE_NOTHROW(ker.finalizeKernel(stat)); }

    SECTION("setCallOnlyOnce prevents multiple executions") {
        ker.setCallOnlyOnce();
        ker.executeKernel(stat);
        ker.executeKernel(stat);  // Should not increment twice
        REQUIRE(ker.count_exec == 1);
        REQUIRE(ker.once_called == true);
    }
}

TEST_CASE("Kernel RuleSet and Metadata", "[Kernel][metadata]") {
    TestKernel ker("MetaKernel");
    auto       rs = createTestRuleSet();

    SECTION("setRuleSet and getRuleSet work correctly") {
        ker.setRuleSet(rs);
        REQUIRE(ker.getRuleSet() == rs);
        REQUIRE(ker._ruleset == rs);
    }

    SECTION("getType returns default 0 for base kernel") {
        REQUIRE(ker.getType() == 0);
        REQUIRE(ker.kernel_type == 0);
    }

    SECTION("operator== compares kernel IDs") {
        TestKernel ker1;
        TestKernel ker2;
        TestKernel ker3(ker1.kernel_name);  // Same name but different ID

        REQUIRE(ker1 == ker1);     // Self comparison
        REQUIRE(!(ker1 == ker2));  // Different IDs
        REQUIRE(!(ker1 == ker3));  // Same name, different ID
    }
}

TEST_CASE("Kernel Information String", "[Kernel][info]") {
    TestKernel   ker("InfoKernel");
    const double test_time  = 42.75;
    const int    test_layer = 2;
    const int    test_depth = 5;
    const int    test_align = 64;

    std::string info = ker.generateInformationString(test_time, test_layer, test_depth, test_align);

    SECTION("Information string contains expected metadata") {
        REQUIRE(info.find("InfoKernel") != std::string::npos);
        REQUIRE(info.find("42.75") != std::string::npos);
        REQUIRE(info.find("2") != std::string::npos);
        REQUIRE(info.find("5") != std::string::npos);
        REQUIRE(info.find("64") != std::string::npos);
    }
}

TEST_CASE("Kernel Timing and Flags", "[Kernel][flags]") {
    TestKernel ker("FlagKernel");

    SECTION("setTiming enables timing flag") {
        ker.setTiming(true);
        REQUIRE(ker.is_timing == true);

        ker.setTiming(false);
        REQUIRE(ker.is_timing == false);
    }

    SECTION("enable_call_child flag controls child execution") {
        REQUIRE(ker.enable_call_child == true);  // Default value

        ker.enable_call_child = false;
        auto child            = std::make_shared<TestKernel>("Child");
        ker.appendChild(child);

        Status stat;
        ker.executeKernel(stat);
        REQUIRE(child->count_exec == 0);  // Child not executed when flag is false
    }
}
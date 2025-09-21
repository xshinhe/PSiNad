#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"
#include "psnd/chem.h"

using namespace chem;
using namespace Catch;

TEST_CASE("Element label retrieval", "[chem][getElemLabel]") {
    SECTION("Valid atomic numbers") {
        REQUIRE(getElemLabel(1) == "H");
        REQUIRE(getElemLabel(8) == "O");
        REQUIRE(getElemLabel(6) == "C");
        REQUIRE(getElemLabel(92) == "U");
        REQUIRE(getElemLabel(118) == "OG");
    }

    SECTION("Invalid atomic numbers") {
        REQUIRE_THROWS_AS(getElemLabel(0), std::runtime_error);
        REQUIRE_THROWS_AS(getElemLabel(119), std::runtime_error);
        REQUIRE_THROWS_AS(getElemLabel(-5), std::runtime_error);
    }
}

TEST_CASE("Element index retrieval", "[chem][getElemIndex]") {
    SECTION("Valid labels (case insensitivity)") {
        REQUIRE(getElemIndex("H") == 1);
        REQUIRE(getElemIndex("h") == 1);
        REQUIRE(getElemIndex("O") == 8);
        REQUIRE(getElemIndex("o") == 8);
        REQUIRE(getElemIndex("C") == 6);
        REQUIRE(getElemIndex("U") == 92);
        REQUIRE(getElemIndex("OG") == 118);
        REQUIRE(getElemIndex("og") == 118);
    }

    SECTION("Invalid labels") {
        REQUIRE(getElemIndex("") == 0);
        REQUIRE(getElemIndex("X") == 0);
        REQUIRE(getElemIndex("ABC") == 0);
        REQUIRE(getElemIndex("123") == 0);
    }
}

TEST_CASE("Element mass retrieval", "[chem][getElemMass]") {
    SECTION("Valid atomic numbers") {
        REQUIRE(getElemMass(1) == Approx(1.007947));
        REQUIRE(getElemMass(8) == Approx(15.99943));
        REQUIRE(getElemMass(6) == Approx(12.01078));
        REQUIRE(getElemMass(92) == Approx(238.028913));
        REQUIRE(getElemMass(118) == Approx(294.0));
    }

    SECTION("Invalid atomic numbers") {
        REQUIRE_THROWS_AS(getElemMass(0), std::runtime_error);
        REQUIRE_THROWS_AS(getElemMass(119), std::runtime_error);
        REQUIRE_THROWS_AS(getElemMass(-3), std::runtime_error);
    }
}

TEST_CASE("ElementList validation", "[chem][ElementList]") {
    SECTION("Check list size") {
        REQUIRE(elem::ElementList.size() == elem::max_number + 1);  // +1 for Z=0 (NU)
    }

    SECTION("Check specific elements") {
        // Hydrogen (Z=1)
        REQUIRE(elem::ElementList[1].label == "H");
        REQUIRE(elem::ElementList[1].Z == 1);
        REQUIRE(elem::ElementList[1].mass == Approx(1.007947));

        // Carbon (Z=6)
        REQUIRE(elem::ElementList[6].label == "C");
        REQUIRE(elem::ElementList[6].Z == 6);
        REQUIRE(elem::ElementList[6].mass == Approx(12.01078));

        // Uranium (Z=92)
        REQUIRE(elem::ElementList[92].label == "U");
        REQUIRE(elem::ElementList[92].Z == 92);
        REQUIRE(elem::ElementList[92].mass == Approx(238.028913));

        // Null element (Z=0)
        REQUIRE(elem::ElementList[0].label == "NU");
        REQUIRE(elem::ElementList[0].Z == 0);
        REQUIRE(elem::ElementList[0].mass == Approx(0.0));
    }
}

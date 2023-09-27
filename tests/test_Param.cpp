#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "../core/Param.h"

using namespace PROJECT_NS;

TEST_CASE("Check Param : fromString", "[Param]") {
    std::string param_str = "{\"x\": 1, \"y\": 2.5}";
    Param P(param_str, Param::fromString);

    REQUIRE(P.get<int>("x") == 1);
    REQUIRE(P.get<double>("y") == 2.5f);
}

TEST_CASE("Check Param : fromFile", "[Param]") {
    std::string param_str = "{\"x\": 1, \"y\": 2.5}";
    std::ofstream ofs("/tmp/tmp.json");
    ofs << param_str;
    ofs.close();

    Param P("/tmp/tmp.json", Param::fromFile);

    REQUIRE(P.get<int>("x") == 1);
    REQUIRE(P.get<double>("y") == 2.5f);
}

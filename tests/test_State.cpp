#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include "../core/DataSet.h"
#include "../generate/version.h"

using namespace kids;

template <typename T>
class Test1_Fixture {
   private:
    int size1 = 1;
    int size2 = 10;

   protected:
    details::Node<T>* pN;

   protected:
    int test_allocator_and_deallocator() {
        try {
            pN = new details::Node<T>(size1);
        } catch (std::runtime_error& e) { return 1; }

        try {
            delete pN;
        } catch (std::runtime_error& e) { return 1; }

        return 0;
    }
};

TEMPLATE_TEST_CASE_METHOD(Test1_Fixture, "Check Node : allocator/deallocator", "[class][template]", int, double,
                          std::complex<double>) {
    REQUIRE(Test1_Fixture<TestType>::test_allocator_and_deallocator() == 0);
}

template <typename T>
class Test2_Fixture {
   private:
    int size1 = 2;
    int size2 = 10;

   protected:
    DataSet* pS;

   protected:
    int test_allocator_and_deallocator() {
        try {
            pS = new DataSet();
        } catch (std::runtime_error& e) { return 1; }

        try {
            delete pS;
        } catch (std::runtime_error& e) { return 1; }
        return 0;
    }

    int test_def() {
        try {
            pS = new DataSet();
            pS->def<T>("data.x.1.1", size1);
            pS->def<T>("data.x.2.2", size2);
            pS->def<int>("int.1.11.1.1.1.1", size1);
            pS->def<double>("double.112113", size1);
            delete pS;
        } catch (std::runtime_error& e) { return 1; }
        return 0;
    }

    int test_undef() {
        try {
            pS = new DataSet();
            pS->def<T>("data.1", size1);
            pS->undef("data.1");
            delete pS;
        } catch (std::runtime_error& e) { return 1; }
        return 0;
    }

    int test_reg() {
        try {
            pS = new DataSet();
            pS->def<T>("data.1.1", size1);
            T* ptr = pS->reg<T>("data.1.1", size1);
            delete pS;
        } catch (std::runtime_error& e) { return 1; }
        return 0;
    }
};

TEMPLATE_TEST_CASE_METHOD(Test2_Fixture, "Check DataSet : allocator/deallocator", "[class][template]", int, double,
                          std::complex<double>) {
    REQUIRE(Test2_Fixture<TestType>::test_allocator_and_deallocator() == 0);
}

TEMPLATE_TEST_CASE_METHOD(Test2_Fixture, "Check DataSet : def()", "[class][template]", int, double,
                          std::complex<double>) {
    REQUIRE(Test2_Fixture<TestType>::test_def() == 0);
}

TEMPLATE_TEST_CASE_METHOD(Test2_Fixture, "Check DataSet : undef()", "[class][template]", int, double,
                          std::complex<double>) {
    REQUIRE(Test2_Fixture<TestType>::test_undef() == 0);
}

TEMPLATE_TEST_CASE_METHOD(Test2_Fixture, "Check DataSet : reg()", "[class][template]", int, double,
                          std::complex<double>) {
    REQUIRE(Test2_Fixture<TestType>::test_reg() == 0);
}

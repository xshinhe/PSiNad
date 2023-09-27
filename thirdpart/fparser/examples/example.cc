// Simple example file for the function parser
// ===========================================

/* When running the program, try for example with these values:

f(x) = x^2
min x: -5
max x: 5
step: 1

*/

#include <complex>
#include <iostream>
#include <string>

#include "../fparser.hh"

int main() {
    using T = std::complex<double>;

    std::string function;
    T minx, maxx, step;

    std::vector<FunctionParserBase<T>> GLOB;

    GLOB.push_back(FunctionParserBase<T>());

    FunctionParserBase<T> fparser;

    fparser.AddConstant("pi", 3.1415926535897932);

    int res = fparser.Parse("x^2 + 1i * x*y", "x,y");

    std::cout << "min x: ";
    std::cin >> minx;
    std::cout << "max x: ";
    std::cin >> maxx;
    std::cout << "step: ";
    std::cin >> step;
    if (std::cin.fail()) return 0;


    std::vector<T> vals = {std::complex(1.0, 0.0), std::complex(0.0, 1.0)};
    std::cout << "f("
              << ") = " << fparser.Eval(vals.data()) << std::endl;

    vals[1] = std::complex(1.0, 0.0);
    std::cout << "f("
              << ") = " << fparser.Eval(vals.data()) << std::endl;

    return 0;
}

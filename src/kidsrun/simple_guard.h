#ifndef SIMPLE_UTILS_H
#define SIMPLE_UTILS_H

#include <type_traits>

class Simple_Guard final {
   public:
    std::size_t istart;
    std::size_t iend;
    std::size_t BEGIN;
    std::size_t TOTAL;

    Simple_Guard(std::size_t BEGIN, std::size_t TOTAL);
    ~Simple_Guard();
};

#endif  // SIMPLE_UTILS_H

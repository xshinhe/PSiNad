#include "simple_guard.h"

Simple_Guard::Simple_Guard(std::size_t BEGIN, std::size_t TOTAL) : BEGIN{BEGIN}, TOTAL{TOTAL} {
    istart = BEGIN;
    iend   = BEGIN + TOTAL;
}

Simple_Guard::~Simple_Guard(){};

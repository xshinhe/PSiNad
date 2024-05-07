#ifndef KIDS_FMT_H
#define KIDS_FMT_H

#include <iomanip>
#include <iostream>

namespace PROJECT_NS {

/**
 * control the io printing format
 */
constexpr inline int FMT_WIDTH(int X) { return X + 7; }
#define FMT(X)                                                            \
    " " << std::setiosflags(std::ios::scientific) /*scientific notation*/ \
        << std::setprecision(X)                   /*precision*/           \
        << std::right                             /*alignment*/           \
        << std::setw(FMT_WIDTH(X))                /*width of text*/

/**
 * show the location information for debug
 */
#define LOC() (std::string(basename(__FILE__)) + ":" + std::to_string(__LINE__))

};  // namespace PROJECT_NS

#endif  // KIDS_FMT_H

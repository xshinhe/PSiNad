#ifndef NS_H
#define NS_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Shape.h"

namespace PROJECT_NS {

namespace DATASET {

extern std::vector<void*> _VAR_LIST;

template <class T>
class VARIABLE {
   public:
    VARIABLE(const std::string& name, Shape* shape, const std::string& doc) : _name{name}, _shape{shape}, _doc{doc} {
        std::string::size_type pos = 0;
        while ((pos = _name.find("::")) != std::string::npos) _name.replace(pos, 2, ".");
        // _VARS.push_back(this);
        _VAR_LIST.push_back((void*) this);
    }

    std::string help() { return _doc; }
    std::string name() { return _name; }
    std::size_t size() { return (*_shape)(); }

   private:
    std::string _name;
    T* _data;
    Shape* _shape;
    std::string _doc;
    bool allocated = false;
};

namespace dimension {
extern VARIABLE<double> x;  ///< No. of MonteCarlo
};                          // namespace dimension

namespace integrator {
namespace baoab {
extern VARIABLE<double> x;
extern VARIABLE<double> p;
extern VARIABLE<double> m;
extern VARIABLE<double> f;
};  // namespace baoab
};  // namespace integrator

};  // namespace DATASET

};  // namespace PROJECT_NS

#endif  // NS_H

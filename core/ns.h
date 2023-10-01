#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "DataSet.h"

namespace PROJECT_NS {

namespace DATASET {

static std::vector<void*> _VAR_LIST;

template <class T>
class VARIABLE {
   public:
    VARIABLE(const std::string& name, const std::string& doc) : _name{name}, _doc{doc} {
        std::string::size_type pos = 0;
        while ((pos = _name.find("::")) != std::string::npos) _name.replace(pos, 2, ".");
        _VAR_LIST.push_back((void*) this);
    }

    T* init_data(DataSet* DS = nullptr, int size = 1) {
        if (!allocated && DS != nullptr) _data = DS->reg<T>(_name, size);
        allocated = true;
        return _data;
    }

    std::string help() { return _doc; }
    std::string name() { return _name; }

   private:
    std::string _name;
    T* _data;
    std::string _doc;
    bool allocated = false;
};

namespace dimension {
extern VARIABLE<int> M;     ///< No. of MonteCarlo
extern VARIABLE<int> P;     ///< No. of MonteCarlo
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

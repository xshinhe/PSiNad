%module(directors="1") kids
%include "factory.i"

%include "std_string.i"
%include "std_iostream.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_set.i"
%include "std_vector.i"

namespace std {
  %template(pairii) pair<int,int>;
  %template(vectord) vector<double>;
  %template(vectorddd) vector< vector< vector<double> > >;
  %template(vectori) vector<int>;
  %template(vectorii) vector < vector<int> >;
  %template(vectorpairii) vector< pair<int,int> >;
  %template(vectorstring) vector<string>;
  %template(mapstringstring) map<string,string>;
  %template(mapstringdouble) map<string,double>;
  %template(mapii) map<int,int>;
  %template(seti) set<int>;
};

%include "typemaps.i"
%include "windows.i"

%{
#define SWIG_FILE_WITH_INIT

#include <sstream>

#include <exception>
#include <fstream>
#include "kids/phys.h"
#include "kids/chem.h"
#include "kids/Param.h"
#include "kids/DataSet.h"
#include "kids/RuleSet.h"
#include "kids/Kernel.h"

using namespace PROJECT_NS;

%}

%feature("autodoc", "0");
%nodefaultctor;

// %include features.i
// %include KIDS_docstring.i
// %include KIDSSwigHeaders.i

%pythoncode %{
  # when we import * from the python module, we only want to import the
  # actual classes, and not the swigregistration methods, which have already
  # been called, and are now unneeded by the user code, and only pollute the
  # namespace
  __all__ = [k for k in locals().keys() if not (k.endswith('_swigregister') or k.startswith('_'))]
%}

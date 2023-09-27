#include "interf_pythonff.h"

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

// #define PyRawRefered(Arr_ptr) Arr_ptr    // it will do copy (safe for data)
#define PyRawRefered(Arr_ptr) \
    Arr_ptr, py::capsule(Arr_ptr, [](void* _void_ptr) { ; })  // zero-copy when pass to py::array

py::module_ pymod_global;
py::object pyfun_parm;
py::object pyfun_init;
py::object pyfun_npes;
py::object pyfun_epes;

PythonFF_ForceField::PythonFF_ForceField(const Param& iparm) : Nad_ForceField(iparm) {
    pymod_file = Param_GetT(std::string, parm, "pymod_file", "py_force.py");

    std::string relative_path = utils::ParseFilePath(pymod_file);
    std::string file_name     = utils::ParseFileName(pymod_file);
    pymod_file                = file_name;

    tag = name() + pymod_file + "_" + tag;

    py::initialize_interpreter();

    auto sys = py::module_::import("sys");
    sys.attr("path").attr("append")(relative_path);

    pymod_global = py::module_::import(pymod_file.c_str());  // only import once

    pyfun_parm = pymod_global.attr("ForceField_parm");
    pyfun_init = pymod_global.attr("ForceField_init");
    pyfun_npes = pymod_global.attr("ForceField_npes");
    pyfun_epes = pymod_global.attr("ForceField_epes");

    pyfun_parm(parm.dump(4, ' '));
};

PythonFF_ForceField::~PythonFF_ForceField() { py::finalize_interpreter(); };

int PythonFF_ForceField::ForceField_init(num_real* nr, num_real* np, num_real* nm, num_complex* erho, num_complex* eeac,
                                         int& eocc, const int& rdim, const int& fdim, const int& itraj) {
    plFunction();
    {
        const py::tuple& res_list                      = pyfun_init(rdim, fdim, itraj);
        const py::array_t<num_real>& nr_array_out      = py::cast<py::array_t<num_real>>(res_list[0]);
        const py::array_t<num_real>& np_array_out      = py::cast<py::array_t<num_real>>(res_list[1]);
        const py::array_t<num_real>& nm_array_out      = py::cast<py::array_t<num_real>>(res_list[2]);
        const py::array_t<num_complex>& erho_array_out = py::cast<py::array_t<num_complex>>(res_list[3]);
        const py::array_t<num_complex>& eeac_array_out = py::cast<py::array_t<num_complex>>(res_list[4]);
        eocc                                           = py::cast<int>(res_list[5]);

        for (int i = 0; i < rdim; ++i) nr[i] = nr_array_out.data()[i];
        for (int i = 0; i < rdim; ++i) np[i] = np_array_out.data()[i];
        for (int i = 0; i < rdim; ++i) nm[i] = nm_array_out.data()[i];
        for (int i = 0; i < fdim * fdim; ++i) erho[i] = erho_array_out.data()[i];
        for (int i = 0; i < fdim; ++i) eeac[i] = eeac_array_out.data()[i];
    }
    return 0;
}

int PythonFF_ForceField::ForceField_npes(num_real* V, num_real* dV, num_real* ddV, num_real* R, num_real* P,
                                         const int& flag, const int& rdim) {
    plFunction();
    {
        // method 1: all memory is accessiable in python code, good for expert but may not safe!
        // the python function can by binded with @jit(nopython=True)
        int NNuse = (flag > 1) ? NN : 0;

        py::array_t<num_real> R_array(rdim, PyRawRefered(R));       // zero-copy, in & out
        py::array_t<num_real> P_array(rdim, PyRawRefered(P));       // zero-copy, in & out
        py::array_t<num_real> V_array(1, PyRawRefered(V));          // zero-copy, in & out
        py::array_t<num_real> dV_array(rdim, PyRawRefered(dV));     // zero-copy, in & out
        py::array_t<num_real> ddV_array(NNuse, PyRawRefered(ddV));  // zero-copy, in & out
        pyfun_npes(V_array, dV_array, ddV_array, R_array, P_array, flag, rdim);
    }
    /*{
        // method 2: return a tuple from python. (only a little loss of efficiency, not so mush)
        py::array_t<num_real> R_array(rdim, PyRawRefered(R));  // zero-copy, in & out
        py::array_t<num_real> P_array(rdim, PyRawRefered(P));  // zero-copy, in & out
        const py::tuple& res_list                  = pyfun_npes(R_array, P_array, flag, rdim);
        const py::array_t<num_real>& V_array_out   = py::cast<py::array_t<num_real>>(res_list[0]);
        const py::array_t<num_real>& dV_array_out  = py::cast<py::array_t<num_real>>(res_list[1]);
        const py::array_t<num_real>& ddV_array_out = py::cast<py::array_t<num_real>>(res_list[2]);
        for (int i = 0; i < 1; ++i) V[i] = V_array_out.data()[i];
        if (flag > 0) {
            for (int i = 0; i < rdim; ++i) dV[i] = dV_array_out.data()[i];
        }
        if (flag > 1) {
            for (int i = 0; i < rdim * rdim; ++i) ddV[i] = ddV_array_out.data()[i];
        }
    }*/
    return 0;
}

int PythonFF_ForceField::ForceField_epes(num_real* V, num_real* dV, num_real* ddV, num_real* R, const int& flag,
                                         const int& rdim, const int& fdim) {
    plFunction();
    {
        // method 1: all memory is accessiable in python code, good for expert but may not safe!
        int NNFFuse = (flag > 1) ? NNFF : 0;
        py::array_t<num_real> R_array(rdim, PyRawRefered(R));         // zero-copy
        py::array_t<num_real> V_array(FF, PyRawRefered(V));           // zero-copy, in & out
        py::array_t<num_real> dV_array(NFF, PyRawRefered(dV));        // zero-copy, in & out
        py::array_t<num_real> ddV_array(NNFFuse, PyRawRefered(ddV));  // zero-copy, in & out
        pyfun_epes(V_array, dV_array, ddV_array, R_array, flag, rdim, fdim);
    }
    /*{
        py::array_t<num_real> R_array(rdim, PyRawRefered(R));  // zero-copy
        const py::tuple& res_list                  = pyfun_epes(R_array, flag, rdim, fdim);
        const py::array_t<num_real>& V_array_out   = py::cast<py::array_t<num_real>>(res_list[0]);
        const py::array_t<num_real>& dV_array_out  = py::cast<py::array_t<num_real>>(res_list[1]);
        const py::array_t<num_real>& ddV_array_out = py::cast<py::array_t<num_real>>(res_list[2]);
        for (int i = 0; i < fdim * fdim; ++i) V[i] = V_array_out.data()[i];
        if (flag > 0) {
            for (int i = 0; i < rdim * fdim * fdim; ++i) dV[i] = dV_array_out.data()[i];
        }
        if (flag > 1) {
            for (int i = 0; i < rdim * rdim * fdim * fdim; ++i) ddV[i] = ddV_array_out.data()[i];
        }
    }*/
    return 0;
}

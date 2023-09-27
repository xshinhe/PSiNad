#ifndef Context_H
#define Context_H

#include "../thirdpart/HighFive/include/highfive/H5DataSet.hpp"
#include "../thirdpart/HighFive/include/highfive/H5DataSpace.hpp"
#include "../thirdpart/HighFive/include/highfive/H5File.hpp"
#include "../thirdpart/HighFive/include/highfive/H5FileDriver.hpp"
#include "../thirdpart/HighFive/include/highfive/H5Group.hpp"
#include "mpi_utils.h"  // parallel io within context

#define DUMP(pfile, data, size) hdump(pfile, #data, data, size)

#define LOAD(pfile, data, size) hload(pfile, #data, data, size)

using namespace HighFive;

template <typename T>
int hdump(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size) {
    if (pfile->exist(name)) {
        // no type & dimensions check! @unsafe
        auto ds = pfile->getDataSet(name);
        ds.select({std::size_t(mpi_rank), 0}, {1, size}).write(data);
    } else {
        auto ds = pfile->createDataSet<T>(name, DataSpace({std::size_t(mpi_nprocs), size}));
        ds.select({std::size_t(mpi_rank), 0}, {1, size}).write(data);
    }
    return 0;
}

template <typename T>
int hload(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size) {
    try {
        auto ds   = pfile->getDataSet(name);
        auto dims = ds.getDimensions();
        CHECK_EQ(dims[0], mpi_nprocs);
        CHECK_EQ(dims[1], size);

        ds.select({std::size_t(mpi_rank), 0}, {1, size}).read(data);
    } catch (HighFive::Exception& e) { LOG(FATAL); }
    return 0;
}
// Using MPI rank to manage data load & dump is too limited and inconvenient for communication
//  between jobs with different MPI procs. Thus I overload this two API  
template <typename T>
int hdump_extend(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size, const int& DIM0, const int& rank) {
    if (pfile->exist(name)) {
        // no type & dimensions check! @unsafe
        auto ds = pfile->getDataSet(name);
        ds.select({std::size_t(rank), 0}, {1, size}).write(data);
    } else {
        auto ds = pfile->createDataSet<T>(name, DataSpace({std::size_t(DIM0), size}));
        ds.select({std::size_t(rank), 0}, {1, size}).write(data);
    }
    return 0;
}

template <typename T>
int hload_extend(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size, const int& DIM0, const int& rank) {
    try {
        auto ds   = pfile->getDataSet(name);
        auto dims = ds.getDimensions();
        CHECK_EQ(dims[0], DIM0);
        CHECK_EQ(dims[1], size);

        ds.select({std::size_t(rank), 0}, {1, size}).read(data);
    } catch (HighFive::Exception& e) { LOG(FATAL); }
    return 0;
}

struct ContextNode {
    enum _type { null_type, int_type, real_type, complex_type };
    union {
        _type _data_type = null_type;
        size_t _data_size;
        int* _data_int;
        num_real* _data_real;
        num_complex* _data_complex;
    } data;

    ContextNode(int* data, const size_t& data_size) {
        _data._data_int  = new int[data_size];
        _data._data_type = int_type;
        for (int i = 0; i < data_size; ++i) _data._data_int[i] = data[i];
    }

    ContextNode(num_real* data, const size_t& data_size) {
        _data._data_real = new num_real[data_size];
        _data._data_type = real_type;
        for (int i = 0; i < data_size; ++i) _data._data_real[i] = data[i];
    }

    ContextNode(num_complex* data, const size_t& data_size) {
        _data._data_complex = new num_complex[data_size];
        _data._data_type    = complex_type;
        for (int i = 0; i < data_size; ++i) _data._data_complex[i] = data[i];
    }

    ~ContextNode() {
        switch (_data._data_type) {
            case int_type:
                delete[] _data._data_int;
                break;
            case real_type:
                delete[] _data._data_real;
                break;
            case complex_type:
                delete[] _data._data_complex;
                break;
        }
        _data._data_type = null_type;
    }

    inline void UpdateNode(int* data, const size_t& data_size) {
        assert(_data._data_type == int_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) _data._data_int[i] = data[i];
    }

    inline void UpdateNode(num_real* data, const size_t& data_size) {
        assert(_data._data_type == real_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) _data._data_real[i] = data[i];
    }

    inline void UpdateNode(num_complex* data, const size_t& data_size) {
        assert(_data._data_type == complex_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) _data._data_complex[i] = data[i];
    }

    inline void ReadNode(int* data, const size_t& data_size) {
        assert(_data._data_type == int_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) data[i] = _data._data_int[i];
    }

    inline void ReadNode(num_real* data, const size_t& data_size) {
        assert(_data._data_type == real_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) data[i] = _data._data_real[i];
    }

    inline void ReadNode(num_complex* data, const size_t& data_size) {
        assert(_data._data_type == complex_type);
        assert(_data._data_size == data_size);
        for (int i = 0; i < data_size; ++i) data[i] = _data._data_complex[i];
    }

   protected:
    data _data;
};

class Context {
    Context();
    ~Context();

    inline bool exist(std::string name) { return _dict.find(name) != _dict.end; }

    template <typename T>
    static int dump(Context* pContext, const std::string& name, T* data, const std::size_t& size) {
        if (pContext->exist(name)) {
            auto m = _dict.find(name);
            m->second.UpdateNode(data, size);
        } else {
            pContext->_dict[name] = ContextNode(data, size);
        }
        return 0;
    }

    template <typename T>
    static int load(Context* pContext, const std::string& name, T* data, const std::size_t& size) {
        if (pContext->exist(name)) {
            auto m = _dict.find(name);
            m->second.ReadNode(data, size);
        } else {
            LOG(FATAL);
        }
        return 0;
    }

   protected:
    std::map<std::string, ContextNode> _dict;
}

#endif  // Context_H

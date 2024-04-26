#ifndef HDF5_UTILS_H
#define HDF5_UTILS_H

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5FileDriver.hpp>
#include <highfive/H5Group.hpp>

#include "mpi_utils.h"

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
int hdump_extend(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size, const int& DIM0,
                 const int& rank) {
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
int hload_extend(HighFive::File* pfile, const std::string& name, T* data, const std::size_t& size, const int& DIM0,
                 const int& rank) {
    try {
        auto ds   = pfile->getDataSet(name);
        auto dims = ds.getDimensions();
        CHECK_EQ(dims[0], DIM0);
        CHECK_EQ(dims[1], size);

        ds.select({std::size_t(rank), 0}, {1, size}).read(data);
    } catch (HighFive::Exception& e) { LOG(FATAL); }
    return 0;
}

#endif  // HDF5_UTILS_H

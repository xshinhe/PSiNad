#ifndef MPI_UTILS_H
#define MPI_UTILS_H
#include "mpi.h"

extern int mpi_rank, mpi_nprocs, mpi_init_tag;
extern bool mpi_isroot;
extern int* mpi_range_array;

// extern MPI_Request* mpi_requests;

int mpi_range(const int& idx1, const int& idx2, const int& nprocs, const int& irank, int& ista, int& iend);

int mpi_utils_init();

int mpi_utils_final();

#endif  // MPI_UTILS_H

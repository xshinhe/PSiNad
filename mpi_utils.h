#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include "mpi.h"

class MPI_Guard final {
   public:
    static int rank;
    static int nprocs;
    static bool isroot;

    MPI_Guard();

    static int range(const int& idx1, const int& idx2, int& ista, int& iend);

    ~MPI_Guard();
};



#endif  // MPI_UTILS_H

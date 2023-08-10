#include "mpi_utils.h"

MPI_Guard::MPI_Guard() {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    isroot = (rank == 0);
}

int MPI_Guard::range(const int& idx1, const int& idx2, int& ista, int& iend) {
    int num = (idx2 - idx1) / nprocs;
    ista    = idx1 + rank * num;
    iend    = (rank == nprocs - 1) ? idx2 : ista + num;
    return 0;
}

MPI_Guard::~MPI_Guard() { MPI_Finalize(); }


int MPI_Guard::rank    = 0;
int MPI_Guard::nprocs  = 1;
bool MPI_Guard::isroot = true;

#include "mpi_utils.h"

#include "io_utils.h"

int mpi_rank = 0, mpi_nprocs = 1;

bool mpi_isroot  = true;
int mpi_init_tag = 0;
int* mpi_range_array;

int mpi_range(const int& idx1, const int& idx2, const int& nprocs, const int& irank, int& ista, int& iend) {
    int num = (idx2 - idx1) / nprocs;
    ista    = idx1 + irank * num;
    iend    = (irank == nprocs - 1) ? idx2 : ista + num;
    return 0;
}

int mpi_utils_init() {
    MPI_Init(global::p_argc, global::p_argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_range_array = new int[mpi_nprocs];

    mpi_isroot   = !mpi_rank;
    mpi_init_tag = 0;

    return 0;
}

int mpi_utils_final() {
    // MPI Finalize
    MPI_Finalize();
    delete[] mpi_range_array;
    return 0;
}


class MPI_Guard final {
   public:
    static int rank;
    static int nprocs;
    static bool isroot;

    MPI_Guard() {
        MPI_Init(global::p_argc, global::p_argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        isroot = (rank == 0);
    }

    static int range(const int& idx1, const int& idx2, int& ista, int& iend) {
        int num = (idx2 - idx1) / nprocs;
        ista    = idx1 + rank * num;
        iend    = (rank == nprocs - 1) ? idx2 : ista + num;
        return 0;
    }

    ~MPI_Guard() { MPI_Finalize(); }
};

int MPI_Guard::rank    = 0;
int MPI_Guard::nprocs  = 1;
bool MPI_Guard::isroot = true;

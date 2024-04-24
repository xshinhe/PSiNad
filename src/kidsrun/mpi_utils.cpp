#include "mpi_utils.h"

namespace PROJECT_NS {
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

int MPI_Guard::reduce(const std::vector<std::tuple<kids_dtype, void*, void*, std::size_t>>& info_list) {
    kids_dtype  dtype;
    void *      from_data, *to_data;
    std::size_t ndata;
    for (auto&& info : info_list) {
        std::tie(dtype, from_data, to_data, ndata) = info;
        switch (dtype) {
            case kids_int_type: {
                MPI_Reduce((kids_int*) from_data, (kids_int*) to_data,  //
                           ndata, MPI::INT, MPI_SUM, 0, MPI_COMM_WORLD);
                break;
            }
            case kids_real_type: {
                MPI_Reduce((kids_real*) from_data, (kids_real*) to_data,  //
                           ndata, MPI::DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
                break;
            }
            case kids_complex_type: {
                MPI_Reduce((kids_complex*) from_data, (kids_complex*) to_data,  //
                           ndata, MPI::DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
                break;
            }
        }
    }
    return 0;
}

MPI_Guard::~MPI_Guard() { MPI_Finalize(); }


int  MPI_Guard::rank   = 0;
int  MPI_Guard::nprocs = 1;
bool MPI_Guard::isroot = true;

};  // namespace PROJECT_NS
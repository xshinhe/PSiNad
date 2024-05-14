#include "mpi_utils.h"

#include <tuple>

namespace PROJECT_NS {
MPI_Guard::MPI_Guard(std::size_t TOTAL) : TOTAL{TOTAL} {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Guard::range(0, TOTAL, istart, iend);
    isroot = (rank == 0);
}

int MPI_Guard::range(const size_t& idx1, const size_t& idx2, size_t& ista, size_t& iend) {
    int num = (idx2 - idx1) / nprocs;
    ista    = idx1 + rank * num;
    iend    = (rank == nprocs - 1) ? idx2 : ista + num;
    return 0;
}

/*
void collect(std::vector<std::tuple<void*, void*, kids_dtype, kids_option, size_t, size_t>>& list, size_t ncalls) {
    for (auto& term : list) {
        auto [fromdata_raw, todata_raw, dataType, ruleOption, startidx, datasize] = term;  //
        switch (dataType) {
            case kids_real_type: {
                kids_real* fromdata = (kids_real*) fromdata_raw;
                kids_real* todata   = (kids_real*) todata_raw;
                switch (ruleOption) {
                    case kids_copy: {
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i) todata[i] = fromdata[i];
                        break;
                    }
                    case kids_sum: {
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i) todata[i] += fromdata[i];
                        break;
                    }
                    case kids_mean: {
                        kids_real k1 = ncalls / (kids_real)(ncalls + 1);
                        kids_real k2 = 1.0e0 - k1;
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i)
                            todata[i] = k1 * todata[i] + k2 * fromdata[i];
                        break;
                    }
                }
                break;
            }
            case kids_complex_type: {
                kids_complex* fromdata = (kids_real*) fromdata_raw;
                kids_complex* todata   = (kids_real*) todata_raw;
                switch (ruleOption) {
                    case kids_copy: {
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i) todata[i] = fromdata[i];
                        break;
                    }
                    case kids_sum: {
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i) todata[i] += fromdata[i];
                        break;
                    }
                    case kids_mean: {
                        kids_real k1 = ncalls / (kids_real)(ncalls + 1);
                        kids_real k2 = 1.0e0 - k1;
                        for (int i = startidx, cnt = 0; cnt < datasize; ++cnt, ++i)
                            todata[i] = k1 * todata[i] + k2 * fromdata[i];
                        break;
                    }
                }
                break;
            }
        }
    }
}
*/

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

int MPI_Guard::reduce(const std::tuple<kids_dtype, void*, void*, std::size_t>& info) {
    // kids_dtype  dtype;
    // void *      from_data, *to_data;
    // std::size_t ndata;
    // std::tie(dtype, from_data, to_data, ndata) = info;
    auto [dtype, from_data, to_data, ndata] = info;
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

    return 0;
}


MPI_Guard::~MPI_Guard() { MPI_Finalize(); }


int  MPI_Guard::rank   = 0;
int  MPI_Guard::nprocs = 1;
bool MPI_Guard::isroot = true;

};  // namespace PROJECT_NS
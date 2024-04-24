#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <vector>

#include "kids/Types.h"
#include "mpi.h"

namespace PROJECT_NS {
class MPI_Guard final {
   public:
    static int  rank;
    static int  nprocs;
    static bool isroot;

    MPI_Guard();

    static int range(const int& idx1, const int& idx2, int& ista, int& iend);

    static int reduce(const std::vector<std::tuple<kids_dtype, void*, void*, std::size_t>>& info_list);

    ~MPI_Guard();
};
};  // namespace PROJECT_NS

#endif  // MPI_UTILS_H

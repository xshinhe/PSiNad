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

    std::size_t istart;
    std::size_t iend;
    std::size_t BEGIN;
    std::size_t TOTAL;

    MPI_Guard(std::size_t BEGIN, std::size_t TOTAL);
    ~MPI_Guard();

    static int reduce(const std::tuple<kids_dtype, void*, void*, std::size_t>& info);
    static int reduce(const std::vector<std::tuple<kids_dtype, void*, void*, std::size_t>>& info_list);

   private:
    static int range(const size_t& idx1, const size_t& idx2, size_t& ista, size_t& iend);
};
};  // namespace PROJECT_NS

#endif  // MPI_UTILS_H

#ifndef KIDS_Status_H
#define KIDS_Status_H

namespace PROJECT_NS {

struct Status {
    Status() : succ{true}, stage{0}, mpi_rank{0}, icalc{0} {};
    Status(bool succ, int stage = 0, int mpi_rank = 0, int icalc = 0)
        : succ{succ}, stage{stage}, mpi_rank{mpi_rank}, icalc{icalc} {};
    bool succ     = true;
    int  stage    = 0;
    int  mpi_rank = 0;
    int  icalc    = 0;
};

};  // namespace PROJECT_NS

#endif  // KIDS_Status_H

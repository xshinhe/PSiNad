#ifndef opendf_control_H
#define opendf_control_H

#include "../utils/definitions.h"
#include "../utils/hdf5_utils.h"
#include "solver.h"

class Control {
   public:
    Control(Solver* pS);
    virtual int run() = 0;
    virtual ~Control();

   protected:
    State* pState;
    Statistic* pStatistic;
    Solver* pSolver;
};

class SingleControl : public Control {
   public:
    SingleControl(Solver* pS);
    virtual ~SingleControl();
    virtual int run();
};

class MPIControl : public Control {
   public:
    MPIControl(Solver* pS);
    virtual int run();
};


class MonteCarloControl : public Control {
   public:
    MonteCarloControl(Solver* pS);
    virtual ~MonteCarloControl();
    virtual int run();
};

#endif  // opendf_control_H
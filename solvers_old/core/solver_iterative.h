#ifndef SOLVER_ITERATIVE_H
#define SOLVER_ITERATIVE_H

template <typename>

class IterativeState {
    typedef std::pair<int, char> Pair;

    std::vector<Pair> pairs;

    pairs.push_back(Pair(0, 'c'));
    pairs.push_back(Pair(1, 'a'));
    pairs.push_back(Pair(42, 'b'));
}

template <typename StateTpl>
class IterativeSolver : Solver {
    IterativeSolver();

    virtual int put_state(std::shared_ptr<IterativeState> state);
    virtual std::shared_ptr<IterativeState> get_state();
    virtual int exec_kernel_impl(int stat = -1);

    Solver* serial_child;
};

// #define OPENDF_STATE_DEFINITION(m, )

#endif  // SOLVER_ITERATIVE_H

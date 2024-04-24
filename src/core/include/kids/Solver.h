#ifndef KIDS_Solver_H
#define KIDS_Solver_H

namespace PROJECT_NS {

class Solver {
    std::shared_ptr<Platform> getPlatfrom() { return _platform; }
    std::shared_ptr<Param>    getParam() { return _param; }
    std::shared_ptr<DataSet>  getDataSet() { return _dataset; }
    std::shared_ptr<System>   getSystem() { return _system; }
    std::shared_ptr<Solver>   getSolver() { return _solver; }

   private:
    std::shared_ptr<Platform> _platform;
    std::shared_ptr<Param>    _param;
    std::shared_ptr<DataSet>  _dataset;
    std::shared_ptr<System>   _system;
    std::shared_ptr<Solver>   _solver;
};
};  // namespace PROJECT_NS

#endif  // KIDS_Solver_H
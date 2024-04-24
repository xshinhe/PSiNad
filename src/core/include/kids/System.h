

#ifndef KIDS_System_H
#define KIDS_System_H

namespace PROJECT_NS {

class System {
    std::shared_ptr<Platform> getPlatfrom() { return _platform; }
    std::shared_ptr<Param>    getParam() { return _param; }
    std::shared_ptr<DataSet>  getDataSet() { return _dataset; }
    std::shared_ptr<System>   getSystem() { return _system; }
    std::shared_ptr<Solver>   getSolver() { return _solver; }

   private:
    double*                                      boxvectors;
    double*                                      masses;
    std::vector<std::shared_ptr<ConstraintInfo>> constraints;
    std::vector<std::shared_ptr<Force>>          forces;
    std::vector<std::shared_ptr<VirtualSite>>    virtualSites;
};
};  // namespace PROJECT_NS

#endif  // KIDS_System_H
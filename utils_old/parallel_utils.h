#ifndef OPENDF_PARALLEL_H
#define OPENDF_PARALLEL_H


namespace ParaPolicy {
enum _enum { no_para, calc_para, traj_para };
const std::map<std::string, _enum> _dict = {{"no_para", no_para}, {"calc_para", calc_para}, {"traj_para", traj_para}};
};  // namespace ParaPolicy

class ParaInfo {
   public:
    std::string flag;
    std::string flag;

    ParaInfo(const std::string& flag) { type = ; }
};

#endif OPENDF_PARALLEL_H

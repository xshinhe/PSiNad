#ifndef SOLVER_H
#define SOLVER_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../models/model.h"
#include "../utils/definitions.h"
#include "../utils/hdf5_utils.h"


namespace ParaPolicy {
enum _enum { no_para, calc_para, traj_para };
const std::map<std::string, _enum> _dict = {{"no_para", no_para}, {"calc_para", calc_para}, {"traj_para", traj_para}};
};  // namespace ParaPolicy



class Solver : public Model {
   public:
    Solver(const Param& iparm, Model* pM);

    virtual ~Solver();

    int run();

    /**
     * @brief initialize simulation context (reading & sampling, but no allocation!)
     */
    virtual int init(const int& flag) { return 0; }

    /**
     * @brief finalizing simulation context
     */
    virtual int final(const int& flag) { return 0; }

    /**
     * @brief backup & flush simulation context
     */
    virtual int cache(const int& flag) { return 0; }

   protected:
    std::string save;
    Model* pModel;
    File* pContext;
    int para_type;
};

#endif  // SOLVER_H

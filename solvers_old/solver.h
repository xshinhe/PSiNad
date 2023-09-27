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
    Solver(const std::string& iparm_str, Model* pM) : Solver(Param::parse(iparm_str), pM){};

    virtual ~Solver();

    /**
     * @brief simulation wrapper
     * @details it doesn't provide realization of a method, its child should
     *   override run_impl() & run_parallel()
     * @param argc number of arguments
     * @param argv array of arguments
     * @return status
     */
    int run();

    /**
     * @brief implement for different simulation (non-parallel version)
     * @details overrided in child class
     * @param argc number of arguments
     * @param argv array of arguments
     * @return status
     */
    virtual int run_impl();

    /**
     * @brief implement for different simulation (parallel version)
     * @details overrided in child class
     * @param argc number of arguments
     * @param argv array of arguments
     * @return status
     */
    virtual int run_parallel();

    /**
     * @brief initialize simulation context (reading & sampling, but no allocation!)
     */
    virtual int init(int flag) { return 0; }

    /**
     * @brief finalizing simulation context
     */
    virtual int final(int flag) { return 0; }

    /**
     * @brief backup & flush simulation context
     */
    virtual int cache(int flag) { return 0; }

    std::string save;
    std::string para_flag;

   protected:
    Model* pModel;
    File* pContext_in, *pContext_out, *pContext;
    int para_type;
};

#endif  // SOLVER_H

/**
 * @file basic_Kernels.h
 * @author xshinhe
 * @date 2023-04
 * @brief this file defines basic kernels:
 *      1) for load/dump state;
 *      2) set up random engine;
 */

#ifndef Kernel_DataSetHandles_H
#define Kernel_DataSetHandles_H

#include "../core/Formula.h"
#include "../core/Kernel.h"

namespace PROJECT_NS {

/**
 * @brief Kernel_Load_DataSet load previous data in state stucture if restart
 *  option is specified
 *  It should be added to the first beginning of every solver's builder.
 */
class Kernel_Load_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Load_DataSet"; }

   private:
    DataSet* pDS;    ///< pointer to DataSet
    std::string fn;  ///< filename

    virtual void read_param_impl(Param* P);

    virtual void init_data_impl(DataSet* DS);

    virtual int exec_kernel_impl(int stat = -1);
};

/**
 * @brief Kernel_Dump_DataSet dump current data in state stucture
 *  option is specified
 *  It should be added to the last end (as well as middle) of every solver's
 *  builder
 */
class Kernel_Dump_DataSet : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Dump_DataSet"; }

   private:
    DataSet* pDS;    ///< pointer to DataSet
    std::string fn;  ///< filename

    virtual void read_param_impl(Param* P);

    virtual void init_data_impl(DataSet* S);

    virtual int exec_kernel_impl(int stat = -1);
};

struct Result {
   public:
    int size;
    int frame;
    std::vector<std::string> header;
    std::vector<int> stat;
    std::vector<double> data;
    std::ofstream ofs;

    Result(){};

    Result(Result& res) {
        size   = res.size;
        frame  = res.frame;
        header = res.header;
        stat.resize(frame);
        data.resize(frame * size);
        memset(stat.data(), 0, frame * sizeof(int));
        memset(data.data(), 0, frame * size * sizeof(double));
    }

    void save(const std::string& fname, double t0, double dt, bool with_header = false);
    virtual ~Result();
};

class Kernel_Record final : public Kernel {
   public:
    inline virtual const std::string name() { return "Kernel_Record"; }

    Kernel_Record(){};

    virtual ~Kernel_Record();

    inline static Result& get_sampling() { return _sampling; }
    inline static Result& get_correlation() { return _correlation; }

   private:
    static Result _sampling;
    static Result _correlation;

    std::ofstream ofs_samp;
    std::ofstream ofs_corr;

    int* istep_ptr;
    int* sstep_ptr;
    int* isamp_ptr;
    int* nsamp_ptr;

    std::vector<int> Sampling_ID;  // 1 point sampling
    std::vector<std::string> Sampling_STR;

    std::vector<int> Correlation_ID1;  // 2 points correlation: zero time
    std::vector<int> Correlation_ID2;  // 2 points correlation: t time
    std::vector<std::string> Correlation_STR1;
    std::vector<std::string> Correlation_STR2;

    bool trace;
    double dt;

    // Result sampling;
    // Result correlation;

    virtual void read_param_impl(Param* P);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_DataSetHandles_H

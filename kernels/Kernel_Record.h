#ifndef Kernel_Record_H
#define Kernel_Record_H

#include "../core/Kernel.h"

namespace PROJECT_NS {

struct Result {
   public:
    double t0, dt;
    int size;
    int frame;
    std::vector<std::string> header;
    std::vector<int> stat;
    std::vector<double> data;
    std::ofstream ofs;

    Result(){};

    Result(Result& res) {
        t0     = res.t0;
        dt     = res.dt;
        size   = res.size;
        frame  = res.frame;
        header = res.header;
        stat.resize(frame);
        data.resize(frame * size);
        memset(stat.data(), 0, frame * sizeof(int));
        memset(data.data(), 0, frame * size * sizeof(double));
    }

    void save(const std::string& fname, int ibegin, int length = -1, bool with_header = false);
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
    double t0, dt, time_unit;
    std::string directory;

    virtual void read_param_impl(Param* P);

    virtual void init_data_impl(DataSet* DS);

    virtual void init_calc_impl(int stat = -1);

    virtual int exec_kernel_impl(int stat = -1);
};

};  // namespace PROJECT_NS


#endif  // Kernel_Record_H

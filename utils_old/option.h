#ifndef OPTION_H
#define OPTION_H

struct TimeSpace {
    double tbegin, tend, dt, smalldt, halfdt, largedt;
    int nstep, sstep, nsamp, nrespa;
}

struct Option {
    int type;
    std::string flag;
};

#define InitialOption(OP, DICT, PM, DEF)                          \
    ({                                                            \
        OP.flag = Param_GetT(std::string, PM, #OP##"_flag", DEF); \
        OP.type = DICT.at(OP.flag);                               \
    })

#endif  // OPTION_H

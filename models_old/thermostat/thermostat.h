#ifndef Thermostat_H
#define Thermostat_H
#include <map>

#include "../model.h"

const int size_sy              = 7;
const num_real wgt_sy[size_sy] = {0.784513610477560f, 0.235573213359357f, -1.17767998417887f, 1.3151863206839063f,
                                  -1.17767998417887f, 0.235573213359357f, 0.784513610477560f};

namespace ThermostatPolicy {
enum _enum { None, Langevin, Andersen, NHC };
const std::map<std::string, _enum> _dict = {{"#none", None}, {"#lang", Langevin}, {"#and", Andersen}, {"#nhc", NHC}};
};  // namespace ThermostatPolicy

class Thermostat : public Model {
   public:
    Thermostat(const Param& iparm);
    Thermostat(const std::string& iparm_str) : Thermostat(Param::parse(iparm_str)){};

    virtual ~Thermostat();

    int init_alloc(const int& N);

    bool dothermo(const int& istep);

    inline void set_gammal(num_real gammal_in) { gammal = gammal_in; }

    virtual int evolve(num_real* nr, num_real* np, num_real* nm, const num_real& dt, const int& N, const int& start = 0,
                       num_real gammal_in = -1.0f);

    int run_NHC(num_real* nr, num_real* np, num_real* nm, const num_real& dt, const int& N, const int& start = 0,
                num_real gammal_in = -1.0f);

    num_real beta;

   protected:
    std::string thermo_flag;
    int thermo_type;

    bool has_allocated = false;

    num_real gammal;

    int ndofs, nchain, nrespa, nthstp;  // nrespa is for NHC

    DEFINE_POINTER_PROTECTED(num_real, nhc_r);
    DEFINE_POINTER_PROTECTED(num_real, nhc_p);
    DEFINE_POINTER_PROTECTED(num_real, nhc_G);  //? unsued @TODO
    DEFINE_POINTER_PROTECTED(num_real, nhc_Q);
};


#endif  // Thermostat_H

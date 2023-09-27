#include "commonflags.h"

DEFINE_string(p, "param.json", "paramemter inputs");
DEFINE_string(s, "ENER,TRAJ,ETRAJ,CORR", "outputers; seperated by comma");
DEFINE_string(r, "", "restart from hdf5 records. left blank to disable it");
DEFINE_string(para_flag, "no_para", "using mpi scheme that each core runs a simulation");

DEFINE_bool(read_seed, false, "keep same random seed from hdf5 records file");

DEFINE_int32(nsave_mpi, 5, "save backups during monte carlo");
DEFINE_int32(nsave_time, 5, "save backups during a simulation");

DEFINE_double(everysave, 1.0f, "interval for saving a context (unit in hours)");

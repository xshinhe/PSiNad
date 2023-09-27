#ifndef COMMONFLAGS_H
#define COMMONFLAGS_H

#include <gflags/gflags.h>
#include <glog/logging.h>

DECLARE_string(p);
DECLARE_string(s);
DECLARE_string(r);
DECLARE_string(para_flag);

DECLARE_bool(read_seed);
DECLARE_bool(open_ofs);

DECLARE_int32(nsave_mpi);
DECLARE_int32(nsave_time);

DECLARE_double(everysave);

#endif  // COMMONFLAGS_H

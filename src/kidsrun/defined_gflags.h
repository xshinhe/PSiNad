#ifndef defined_gflags_H
#define defined_gflags_H

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "kids/Param.h"

using namespace PROJECT_NS;

DECLARE_bool(w);
DECLARE_string(p);
DECLARE_string(handler);
DECLARE_string(d);
DECLARE_string(load);
DECLARE_string(dump);
DECLARE_double(backup_time);
DECLARE_bool(restart);
DECLARE_bool(timing);
DECLARE_bool(profiling);
DECLARE_int32(seed);
DECLARE_int32(BGIDX);

void check_and_sync_from_gflags(std::shared_ptr<Param> PM);

#endif  // defined_gflags_H
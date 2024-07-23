#ifndef defined_gflags_H
#define defined_gflags_H

DECLARE_bool(w);
DECLARE_string(p);
DECLARE_string(handler);
DECLARE_string(d);
DECLARE_string(load);
DECLARE_string(dump);
DECLARE_double(backup_time);
DECLARE_bool(timing);
DECLARE_bool(profiling);

void check_and_sync_from_gflags(std::shared_ptr<Param> PM);

#endif  // defined_gflags_H
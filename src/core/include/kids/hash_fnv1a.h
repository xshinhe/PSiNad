#ifndef KIDS_HASH_FNV1A_H
#define KIDS_HASH_FNV1A_H

#include <cstdint>

namespace utils {

constexpr uint32_t fnv1a(const char* str, uint32_t hash = 2166136261u) {
    return (*str == '\0') ? hash : fnv1a(str + 1, (hash ^ static_cast<uint32_t>(*str)) * 16777619u);
}

constexpr uint32_t hash(const char* str) { return fnv1a(str); }

};  // namespace utils


#endif  // KIDS_HASH_FNV1A_H

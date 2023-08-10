#ifndef OPENDF_EXCEPTION_H
#define OPENDF_EXCEPTION_H


#include "concat.h"

namespace PROJECT_NS {
struct basic_error : public std::runtime_error {
    basic_error(std::string const text) : std::runtime_error(text) {}
};

struct state_undefined_key_error : public basic_error {
    state_undefined_key_error(std::string const text) : basic_error(utils::concat("undefined key error of : ", text)) {}
};

struct state_conflicted_key_error : public basic_error {
    state_conflicted_key_error(std::string const text)
        : basic_error(utils::concat("conflicted key error of : ", text)) {}
};

struct state_mismatched_size_error : public basic_error {
    state_mismatched_size_error(std::string const text)
        : basic_error(utils::concat("mismatched size error of : ", text)) {}
};

struct state_mismatched_type_error : public basic_error {
    state_mismatched_type_error(std::string const text)
        : basic_error(utils::concat("mismatched type error of : ", text)) {}
};

struct state_zero_size_error : public basic_error {
    state_zero_size_error(std::string const text) : basic_error(utils::concat("zero size error of : ", text)) {}
};

struct state_load_error : public basic_error {
    state_load_error(std::string const text) : basic_error(utils::concat("state load error of : ", text)) {}
};

struct state_dump_error : public basic_error {
    state_dump_error(std::string const text) : basic_error(utils::concat("state dump error of : ", text)) {}
};

struct param_mismatched_type_error : public basic_error {
    param_mismatched_type_error(std::string const text)
        : basic_error(utils::concat("mismatched type error of : ", text)) {}
};

struct param_illegal_key_error : public basic_error {
    param_illegal_key_error(std::string const text) : basic_error(utils::concat("illegal key error of : ", text)) {}
};

struct param_warning : public basic_error {
    param_warning(std::string const text) : basic_error(utils::concat("parameter warning of : ", text)) {}
};

struct kernel_execute_sequence_error : public basic_error {
    kernel_execute_sequence_error(std::string const text)
        : basic_error(utils::concat("kernel execute sequence error of : ", text)) {}
};

};  // namespace PROJECT_NS

#endif  // OPENDF_EXCEPTION_H

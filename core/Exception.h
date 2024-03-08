#ifndef OPENDF_EXCEPTION_H
#define OPENDF_EXCEPTION_H


#include "concat.h"

namespace PROJECT_NS {

const std::string MSG_PARAM_BAD_KEY    = "bad key error : ";
const std::string MSG_PARAM_BAD_TYPE   = "bad type error : ";
const std::string MSG_DATASET_BAD_KEY  = "bad key error : ";
const std::string MSG_DATASET_BAD_SIZE = "bad size error : ";
const std::string MSG_DATASET_BAD_TYPE = "bad type error : ";
const std::string MSG_DATASET_ERR_LOAD = "load error : ";
const std::string MSG_DATASET_ERR_DUMP = "dump error : ";

struct basic_error : public std::runtime_error {
    basic_error(std::string const text) : std::runtime_error(text) {}
};

struct param_warning : public basic_error {
    param_warning(std::string const text) : basic_error(utils::concat("parameter warning of : ", text)) {}
};

};  // namespace PROJECT_NS

#endif  // OPENDF_EXCEPTION_H

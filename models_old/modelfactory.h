#ifndef MODELFACTORY_H
#define MODELFACTORY_H

#include "../utils/param_wrapper.h"
#include "model.h"

Model* init_model(const std::string& model_name, const Param& parm);

#endif  // MODELFACTORY_H
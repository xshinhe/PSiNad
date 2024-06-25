#include "kids/Model.h"

namespace PROJECT_NS {

Model::Model(const std::string& customized_name) : Kernel(customized_name){};

void Model::setInputParam_impl(std::shared_ptr<Param> PM){};

void Model::setInputDataSet_impl(std::shared_ptr<DataSet> DS){};

Status& Model::initializeKernel_impl(Status& stat) { return stat; }

Status& Model::executeKernel_impl(Status& stat) { return stat; }

Status& Model::finalizeKernel_impl(Status& stat) { return stat; }

};  // namespace PROJECT_NS

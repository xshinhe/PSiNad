#include "kids/System.h"

namespace PROJECT_NS {

System::System(std::shared_ptr<Model> model, std::shared_ptr<Param> PM, std::shared_ptr<DataSet> DS)
    : _model{model}, _param{PM}, _dataset{DS} {
    // _model->setInputParam(_param);
    // _model->setInputDataSet(_dataset);
    // _dim       = _param->get<int>("D", LOC(), 1);
    // _natom     = _param->get<int>("Natom", LOC(), 1);
    // _nmole     = _param->get<int>("Nmole", LOC(), 1);  // unnecessary
    // int _n_all = _param->get<int>("N");
    // int _n_virt = 0;
    // assert(_n_all == _dim * _natom + _n_virt);
}

};  // namespace PROJECT_NS

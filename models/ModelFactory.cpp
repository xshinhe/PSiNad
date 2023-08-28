#include "ModelFactory.h"

#include "Model_Interf_MNDO.h"
#include "Model_LVCM.h"
#include "Model_MD1D.h"
#include "Model_NAD1D.h"
#include "Model_SystemBath.h"

namespace PROJECT_NS {

std::shared_ptr<Kernel> ModelFactory(const std::string& name) {
    if (false) {
    } else if (name == "HO") {
        return std::shared_ptr<Model_HO>(new Model_HO());
    } else if (name == "SystemBath") {
        return std::shared_ptr<Model_SystemBath>(new Model_SystemBath());
    } else if (name == "LVCM") {
        return std::shared_ptr<Model_LVCM>(new Model_LVCM());
    } else if (name == "NAD1D") {
        return std::shared_ptr<Model_NAD1D>(new Model_NAD1D());
    } else if (name == "Interf_MNDO") {
        return std::shared_ptr<Model_Interf_MNDO>(new Model_Interf_MNDO());
    } else {
        throw std::runtime_error("unknown Model name");
    }
    return nullptr;
}

};  // namespace PROJECT_NS

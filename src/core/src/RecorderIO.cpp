#include "kids/RecorderIO.h"

namespace PROJECT_NS {

std::vector<std::shared_ptr<RecorderIO>>& RecorderIO::getRecorderIOs() {
    static std::vector<std::shared_ptr<RecorderIO>> _GLOBAL;
    return _GLOBAL;
}

};  // namespace PROJECT_NS
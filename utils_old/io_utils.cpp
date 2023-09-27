#include "io_utils.h"

#include <fstream>
#include <map>
#include <sstream>

#include "profiling_utils.h"

namespace global {

int* p_argc;
char*** p_argv;

namespace OFS {
std::map<_enum, bool> _isopen;
};  // namespace OFS

int parse_ostream(const std::string& str) {
    std::istringstream isstr(str);
    char delim = ',';
    std::string flag;
    while (getline(isstr, flag, delim)) {
        try {
            OFS::_isopen[OFS::_dict.at(flag)] = true;
        } catch (std::out_of_range& e) {
            int cnt = 0;
            std::stringstream sstr("");
            sstr << "only support:\n-s=\n";
            for (auto& it : OFS::_dict) {
                sstr << FMT(0) << it.first;
                if (++cnt % 4 == 0) sstr << "\n";
            }
            LOG(FATAL) << sstr.str() << std::endl << flag << " is invalid\n" << e.what();
        }
    }
    getline(isstr, flag, '\n');  // change delim back to '\n'
    return 0;
}

};  // namespace global

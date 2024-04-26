#ifndef IO_UTILS_H
#define IO_UTILS_H
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

#include "commonflags.h"

#define COLOR_R std::string("\033[31m")  // red color
#define COLOR_G std::string("\033[32m")  // green color
#define COLOR_B std::string("\033[34m")  // blue color
#define COLOR_Y std::string("\033[33m")  // yellow color
#define COLOR_E std::string("\033[0m")   // end close

constexpr int FMT_WIDTH_SIZE(int X) { return X + 9; }
#define FMT(X) " " << std::setiosflags(std::ios::scientific) << std::setprecision(X) << std::setw(FMT_WIDTH_SIZE(X))

namespace utils {

inline int removeFile(std::string& filename) { return remove(filename.c_str()); }

inline void clearFile(std::string& filename) { std::ofstream clear(filename, std::ios::trunc); }

inline void closeOFS(std::ofstream& ofs) {
    if (ofs.is_open()) ofs.close();
}

inline bool isFileExists(const std::string& name) { return std::ifstream{name.c_str()}.good(); }

#ifdef _WIN32
#define DELIM_SEP '\\'
#else
#define DELIM_SEP '/'
#endif

inline std::string ParseFilePath(const std::string& s) {
    std::string::size_type ipos = s.find_last_of(DELIM_SEP);
    return s.substr(0, ipos);
}

inline std::string ParseFileName(const std::string& s) {
    std::string::size_type ipos = s.find_last_of(DELIM_SEP) + 1;
    std::string filename        = s.substr(ipos, s.length() - ipos);
    return filename.substr(0, filename.rfind("."));
}

inline std::string ParseFileExtension(const std::string& s) {
    std::string::size_type ipos = s.find_last_of(DELIM_SEP) + 1;
    std::string filename        = s.substr(ipos, s.length() - ipos);
    std::string::size_type kpos = filename.rfind(".") + 1;
    return filename.substr(kpos, filename.length() - kpos);
}

inline int copyfile_from_to(const std::string& from, const std::string& to) {
    if (from == to) return 0;
    std::ifstream in(from, std::ios_base::in | std::ios_base::binary);
    std::ofstream out(to, std::ios_base::out | std::ios_base::binary);

    const static int BUF_SIZE = 4096;
    char buf[BUF_SIZE];

    do {
        in.read(&buf[0], BUF_SIZE);
        out.write(&buf[0], in.gcount());
    } while (in.gcount() > 0);
    in.close();
    out.close();
    return 0;
}

};  // namespace utils

namespace global {

extern int* p_argc;
extern char*** p_argv;

namespace OFS {
enum _enum { ENER, SAMP, TRAJ, ESAMP, ETRAJ, CORR };
const std::map<std::string, _enum> _dict = {
    {"ENER", ENER}, {"SAMP", SAMP}, {"TRAJ", TRAJ}, {"ESAMP", ESAMP}, {"ETRAJ", ETRAJ}, {"CORR", CORR},
};
extern std::map<_enum, bool> _isopen;
};  // namespace OFS

// extern std::map<OFS_ENUM, bool> OFS_IS_OPEN_MAP;

int parse_ostream(const std::string& str);

inline bool ofs_is_open(const global::OFS::_enum& test) {
    return (global::OFS::_isopen.find(test) != global::OFS::_isopen.end());
}

inline int del_ostream() { return 0; }

};  // namespace global


#endif  // IO_UTILS_H
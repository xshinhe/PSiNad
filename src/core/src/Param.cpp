#include "kids/Param.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#include "kids/Exception.h"
#include "kids/Types.h"

namespace PROJECT_NS {

Param::Param(const std::string &input, LoadOption option) {
    pj = std::shared_ptr<JSON>(new JSON());
    // clang-format off
    switch (option) {
        case fromFile: {
            std::ifstream ifs(input);
            try {
                ifs >> *pj;
            } catch (const JSON_Exception& e) {
                std::cout << "Invalid file: " << input;
            }
            ifs.close();
            break;
        }
        case fromString: {
            std::stringstream sstr(input);
            try {
                sstr >> *pj;
            } catch (const JSON_Exception& e) {
                std::cout << "Invalid string: " << input;
            }
            break;
        }
    }
    // clang-format on
}

bool Param::has_key(const std::string &key) { return !(pj->count(key) == 0); }

std::shared_ptr<configor::json> Param::pjson() { return pj; }

std::string Param::repr() { return pj->dump(4, ' '); }


template <typename T, bool Require = false>
T Param::get(const std::string &key, const std::string &loc, const phys::dimension7 &qdim, const T &default_value) {
    if (!has_key(key)) {
        if (Require) {
            throw kids_error(                              //
                utils::concat(loc,                         //
                              " Type<", as_str<T>(), ">",  //
                              " Key/", key, "/",           //
                              " : Illegal default")        //
            );
            return T();
        } else {
            try {
                throw param_warning(                           //
                    utils::concat(loc,                         //
                                  " Type<", as_str<T>(), ">",  //
                                  " Key/", key, "/",           //
                                  " : Use default ",           //
                                  default_value)               //
                );
            } catch (param_warning &w) { std::cerr << w.what() << "\n"; }
            return default_value;
        }
    }
    // the case find the key
    switch ((*pj)[key].type()) {
        case configor::config_value_type::string: {
            if (std::is_same<T, std::string>::value) {
                return (*pj)[key].get<T>();
            } else if (std::is_same<T, double>::value) {
                // parse unit
                phys::uval uv   = phys::us::parse((*pj)[key].as_string());
                double     qval = phys::au::as(qdim, uv);

                // conversion by stringstream (stupid)
                T                 q;
                std::stringstream ss;
                ss << std::setiosflags(std::ios::scientific)  //
                   << std::setprecision(32) << qval;          //
                ss >> q;
                return q;
            }
            break;
        }
        case configor::config_value_type::boolean: {
            if (std::is_same<T, bool>::value) return (*pj)[key].get<T>();
            break;
        }
        case configor::config_value_type::number_float: {
            T q;
            if (std::is_same<T, double>::value) {
                q = (*pj)[key].as_float();
                return q;
            }
            break;
        }
        case configor::config_value_type::number_integer: {
            T q;
            if (std::is_same<T, double>::value) {
                q = (*pj)[key].as_float();
                return q;
            } else if (std::is_same<T, int>::value) {
                return (*pj)[key].get<T>();
            }
            break;
        }
    }
    // cannot be adapted to existing conversions
    throw kids_error(                                //
        utils::concat(loc,                           //
                      " Type<", as_str<T>(), ">",    //
                      " Key/", key, "/",             //
                      " Data{",                      //
                      (*pj)[key].dump(4, ' '), "}",  //
                      " : Converting fatal")         //
    );
    return T();
}

/// @{
bool Param::get_bool(const std::string &key, const std::string &loc, const bool &default_value) {
    return get<bool, false>(key, loc, phys::none_d, default_value);
}
bool Param::get_bool(const std::string &key, const std::string &loc) {
    return get<bool, true>(key, loc, phys::none_d, bool());
}

int Param::get_int(const std::string &key, const std::string &loc, const int &default_value) {
    return get<int, false>(key, loc, phys::none_d, default_value);
}
int Param::get_int(const std::string &key, const std::string &loc) {
    return get<int, true>(key, loc, phys::none_d, int());
}

std::string Param::get_string(const std::string &key, const std::string &loc, const std::string &default_value) {
    return get<std::string, false>(key, loc, phys::none_d, default_value);
}
std::string Param::get_string(const std::string &key, const std::string &loc) {
    return get<std::string, true>(key, loc, phys::none_d, std::string());
}

double Param::get_double(const std::string &key, const std::string &loc, const phys::dimension7 &qdim,
                         const double &default_value) {
    return get<double, false>(key, loc, qdim, default_value);
}
double Param::get_double(const std::string &key, const std::string &loc, const double &default_value) {
    return get<double, false>(key, loc, phys::none_d, default_value);
}
double Param::get_double(const std::string &key, const std::string &loc) {
    return get<double, true>(key, loc, phys::none_d, double());
}
/// @}

};  // namespace PROJECT_NS

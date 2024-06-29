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
                throw kids_error("Invalid Json File Format");
            }
            ifs.close();
            break;
        }
        case fromString: {
            std::stringstream sstr(input);
            try {
                sstr >> *pj;
            } catch (const JSON_Exception& e) {
                throw kids_error("Invalid Json String Format");
            }
            break;
        }
    }
    // clang-format on
}

static bool has_key_internal(const Param::JSON &json, const std::string &key) {
    auto        ipos = key.find_first_of(".");
    std::string key1 = (ipos == std::string::npos) ? key : key.substr(0, ipos);
    std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
    if (ipos != std::string::npos) return has_key_internal(json[key1], key2);
    return !(json.count(key) == 0);
}

bool Param::has_key(const std::string &key) { return has_key_internal(*pj, key); }

std::shared_ptr<configor::json> Param::pjson() { return pj; }

std::string Param::repr() { return pj->dump(4, ' '); }

/**
 * @brief get parameter
 * @param json : JSON object
 * @param key : key of the parameter
 * @param loc : location tracer (always leaves it `LOC()`)
 * @param qdim : physics dimension for conversion in double in AU unit
 * @param default_value: default value is avialable
 * @return T type
 */
template <typename T, bool Required = false>
static T get_internal(const Param::JSON &json, const std::string &key, const std::string &loc,  //
                      const phys::dimension7 &qdim, const T &default_value) {
    auto        ipos       = key.find_first_of(".");
    std::string key1       = (ipos == std::string::npos) ? key : key.substr(0, ipos);
    bool        key1_exist = !(json.count(key1) == 0);

    if (!key1_exist && Required) {
        throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Required"));
        return T();
    }
    if (!key1_exist && !Required) {
        try {
            throw param_warning(
                utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Default=", default_value));
        } catch (param_warning &w) { std::cerr << w.what() << "\n"; }
        return default_value;
    }
    if (ipos != std::string::npos) {
        std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
        return get_internal<T, Required>(json[key1], key2, loc, qdim, default_value);
    }

    // the case find the key
    switch (json[key].type()) {
        case configor::config_value_type::string: {
            if (std::is_same<T, std::string>::value) {
                return json[key].get<T>();
            } else if (std::is_same<T, double>::value) {
                phys::uval        uv   = phys::us::parse(json[key].as_string());
                double            qval = phys::au::as(qdim, uv);
                T                 q;
                std::stringstream ss;
                ss << std::setiosflags(std::ios::scientific) << std::setprecision(32) << qval;  // !stupid
                ss >> q;
                return q;
            }
            break;
        }
        case configor::config_value_type::boolean: {
            T q;
            if (std::is_same<T, bool>::value) q = json[key].get<T>();
            return q;
            break;
        }
        case configor::config_value_type::number_float: {
            T q;
            if (std::is_same<T, double>::value) q = json[key].as_float();
            return q;
            break;
        }
        case configor::config_value_type::number_integer: {
            T q;
            if (std::is_same<T, double>::value) q = json[key].as_float();
            if (std::is_same<T, int>::value) q = json[key].get<T>();
            return q;
            break;
        }
    }
    // cannot be adapted to existing conversions
    throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ Data{", json[key].dump(4, ' '), "}",
                                   " : Converting fatal"));
    return T();
}

template <typename T, bool Required = false>
T get(const Param::JSON &json, const std::vector<std::string> &keys, const std::string &loc,  //
      const phys::dimension7 &qdim, const T &default_value) {
    for (auto &key : keys) {
        if (has_key_internal(json, key)) return get_internal<T, Required>(json, key, loc, qdim, default_value);
    }
    if (!Required) return default_value;
    throw kids_error(utils::concat(loc, ": Cannot get parameter for! key = ", keys[0]));
    return T();
}

/// @deprecated
/// @{
// bool Param::get_bool(const std::string &key, const std::string &loc, const bool &default_value) {
//     return get_internal<bool, false>(*pj, key, loc, phys::none_d, default_value);
// }
// bool Param::get_bool(const std::string &key, const std::string &loc) {
//     return get_internal<bool, true>(*pj, key, loc, phys::none_d, bool());
// }

// int Param::get_int(const std::string &key, const std::string &loc, const int &default_value) {
//     return get_internal<int, false>(*pj, key, loc, phys::none_d, default_value);
// }
// int Param::get_int(const std::string &key, const std::string &loc) {
//     return get_internal<int, true>(*pj, key, loc, phys::none_d, int());
// }

// std::string Param::get_string(const std::string &key, const std::string &loc, const std::string &default_value) {
//     return get_internal<std::string, false>(*pj, key, loc, phys::none_d, default_value);
// }
// std::string Param::get_string(const std::string &key, const std::string &loc) {
//     return get_internal<std::string, true>(*pj, key, loc, phys::none_d, std::string());
// }

// double Param::get_real(const std::string &key, const std::string &loc, const phys::dimension7 &qdim,
//                          const double &default_value) {
//     return get_internal<double, false>(*pj, key, loc, qdim, default_value);
// }
// double Param::get_real(const std::string &key, const std::string &loc, const double &default_value) {
//     return get_internal<double, false>(*pj, key, loc, phys::none_d, default_value);
// }
// double Param::get_real(const std::string &key, const std::string &loc) {
//     return get_internal<double, true>(*pj, key, loc, phys::none_d, double());
// }
/// @}

/// new implementation
/// @{
kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc,
                          const kids_bool &default_value) {
    return get<kids_bool, false>(*pj, keys, loc, phys::none_d, default_value);
}
kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc) {
    return get<kids_bool, true>(*pj, keys, loc, phys::none_d, kids_bool());
}

kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc, const kids_int &default_value) {
    return get<kids_int, false>(*pj, keys, loc, phys::none_d, default_value);
}
kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc) {
    return get<kids_int, true>(*pj, keys, loc, phys::none_d, kids_int());
}

std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc,
                              const std::string &default_value) {
    return get<std::string, false>(*pj, keys, loc, phys::none_d, default_value);
}
std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc) {
    return get<std::string, true>(*pj, keys, loc, phys::none_d, std::string());
}

kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc, const phys::dimension7 &qdim,
                          const kids_real &default_value) {
    return get<kids_real, false>(*pj, keys, loc, qdim, default_value);
}
kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc,
                          const kids_real &default_value) {
    return get<kids_real, false>(*pj, keys, loc, phys::none_d, default_value);
}
kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc) {
    return get<kids_real, true>(*pj, keys, loc, phys::none_d, kids_real());
}
/// @}

};  // namespace PROJECT_NS

#include "kids/Param.h"

#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

#include "configor/json.hpp"
#include "kids/Exception.h"
#include "kids/Types.h"
#include "toml11/include/toml.hpp"

namespace PROJECT_NS {

// #ifdef PARAM_USE_JSON

class JsonObject : public Param::DataObject {
   public:
    using DataType  = configor::json;
    using Exception = configor::configor_exception;
    DataType data;
};

class TomlObject : public Param::DataObject {
   public:
    using DataType = toml::value;
    DataType data;
};

Param::Param(const std::string &input, LoadOption option) {
    if (option == fromFile) {
        auto ipos   = input.find_first_of(".");
        auto suffix = input.substr(ipos + 1, input.size());
        if (suffix == "json") {
            impl_t = JSON;
            obj    = std::shared_ptr<DataObject>(new JsonObject());
            std::ifstream ifs(input);
            try {
                ifs >> std::dynamic_pointer_cast<JsonObject>(obj)->data;
            } catch (const configor::configor_exception &e) { throw kids_error("Invalid Json File Format"); }
            ifs.close();
        } else if (suffix == "toml") {
            impl_t = TOML;
            obj    = std::shared_ptr<DataObject>(new TomlObject());
            try {
                std::dynamic_pointer_cast<TomlObject>(obj)->data = toml::parse(input);
            } catch (const toml::exception &e) { throw kids_error("Invalid Toml File Format"); }
        } else {
            throw kids_error("unknown suffix of input! only json & toml are supported.");
        }
    } else if (option == fromString) {
        impl_t = (input[0] == '{') ? JSON : TOML;
        if (impl_t == JSON) {
            obj = std::shared_ptr<DataObject>(new JsonObject());
            std::stringstream sstr(input);
            try {
                sstr >> std::dynamic_pointer_cast<JsonObject>(obj)->data;
            } catch (const configor::configor_exception &e) { throw kids_error("Invalid Json String Format"); }
        } else if (impl_t == TOML) {
            obj = std::shared_ptr<DataObject>(new TomlObject());
            try {
                std::dynamic_pointer_cast<TomlObject>(obj)->data = toml::parse_str(input);
            } catch (const toml::exception &e) { throw kids_error("Invalid Toml String Format"); }
        } else {
            throw kids_error("unknown string format.");
        }
    } else {
        throw kids_error("unknown option for Param");
    }
}

namespace internal {
bool is_true(const JsonObject::DataType &data) { return true; }
bool is_bool(const JsonObject::DataType &data) { return data.is_bool(); }
bool is_int(const JsonObject::DataType &data) { return data.is_integer(); }
bool is_real(const JsonObject::DataType &data) { return data.is_float(); }
bool is_string(const JsonObject::DataType &data) { return data.is_string(); }
bool is_array(const JsonObject::DataType &data) { return data.is_array(); }
bool is_object(const JsonObject::DataType &data) { return data.is_object(); }
bool contains(const JsonObject::DataType &data, const std::string &key) {
    if (is_object(data)) return !(data.count(key) == 0);
    if (is_array(data)) {
        int idx = stoi(key);
        return (0 <= idx && idx < data.size());
    }
    return false;
}
JsonObject::DataType &child(JsonObject::DataType &data, const std::string &key) {
    if (is_object(data)) return data[key];
    if (is_array(data)) return data[stoi(key)];
    throw kids_error("bad child");
    return data;
}
template <class T>
inline T as_type(const JsonObject::DataType &data) {
    return data.get<T>();
}
inline bool        as_bool(const JsonObject::DataType &data) { return data.as_bool(); }
inline int         as_int(const JsonObject::DataType &data) { return data.as_integer(); }
inline double      as_real(const JsonObject::DataType &data) { return data.get<double>(); }
inline std::string as_string(const JsonObject::DataType &data) { return data.as_string(); }
inline std::string dump(const JsonObject::DataType &data) { return data.dump(4, ' '); }

bool is_true(const TomlObject::DataType &data) { return true; }
bool is_bool(const TomlObject::DataType &data) { return data.is_boolean(); }
bool is_int(const TomlObject::DataType &data) { return data.is_integer(); }
bool is_real(const TomlObject::DataType &data) { return data.is_floating(); }
bool is_string(const TomlObject::DataType &data) { return data.is_string(); }
bool is_array(const TomlObject::DataType &data) { return data.is_array(); }
bool is_object(const TomlObject::DataType &data) { return data.is_table(); }
bool contains(const TomlObject::DataType &data, const std::string &key) {
    if (is_object(data)) return data.contains(key);
    if (is_array(data)) {
        int idx = stoi(key);
        return (0 <= idx && idx < data.size());
    }
    return false;
}
TomlObject::DataType &child(TomlObject::DataType &data, const std::string &key) {
    if (is_object(data)) return data[key];
    if (is_array(data)) return data[stoi(key)];
    throw kids_error("bad child");
    return data;
}
template <class T>
inline T as_type(const TomlObject::DataType &data) {
    return toml::get<T>(data);
}
inline bool        as_bool(const TomlObject::DataType &data) { return toml::get<bool>(data); }
inline int         as_int(const TomlObject::DataType &data) { return toml::get<int>(data); }
inline double      as_real(const TomlObject::DataType &data) { return toml::get<double>(data); }
inline std::string as_string(const TomlObject::DataType &data) { return toml::get<std::string>(data); }
inline std::string dump(const TomlObject::DataType &data) { return toml::format(data); }
};  // namespace internal


template <class DType>
static bool apply_internal(DType &data, const std::string &key,  //
                           std::function<bool(const DType &data)> func) {
    auto        ipos = key.find_first_of(".");
    std::string key1 = (ipos == std::string::npos) ? key : key.substr(0, ipos);
    std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
    if (ipos != std::string::npos) {
        if (!internal::contains(data, key1)) return false;
        return apply_internal(internal::child(data, key1), key2, func);
    }
    if (!internal::contains(data, key)) return false;
    return func(internal::child(data, key));
}

bool Param::has_key(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_true);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_true);
    return false;
}

bool Param::is_object(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_object);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_object);
    return false;
}

bool Param::is_array(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_array);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_array);
    return false;
}

bool Param::is_bool(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_bool);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_bool);
    return false;
}

bool Param::is_int(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_int);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_int);
    return false;
}

bool Param::is_real(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_real);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_real);
    return false;
}

bool Param::is_string(const std::string &key) {
    if (impl_t == JSON)
        return apply_internal<JsonObject::DataType>(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
                                                    (bool (*)(const JsonObject::DataType &)) internal::is_string);
    if (impl_t == TOML)
        return apply_internal<TomlObject::DataType>(std::dynamic_pointer_cast<TomlObject>(obj)->data, key,
                                                    (bool (*)(const TomlObject::DataType &)) internal::is_string);
    return false;
}

template <class DType, class T>
static void set_internal(DType &data, const std::string &key, const T &val) {
    auto        ipos = key.find_first_of(".");
    std::string key1 = (ipos == std::string::npos) ? key : key.substr(0, ipos);
    std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
    if (ipos != std::string::npos) {
        if (internal::is_array(data)) {
            int idx1 = stoi(key1);
            if (idx1 >= 0 && idx1 < data.size()) return set_internal(data[idx1], key2, val);
        }
        return set_internal(data[key1], key2, val);
    }
    if (internal::is_array(data)) {
        int idx = stoi(key);
        if (idx >= 0 && idx < data.size()) {
            data[idx] = val;
            return;
        }
    }
    data[key] = val;
}

void Param::set_bool(const std::string &key, bool val) {
    if (impl_t == JSON) set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
    if (impl_t == TOML) set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val);
}

void Param::set_int(const std::string &key, int val) {
    if (impl_t == JSON) set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
    if (impl_t == TOML) set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val);
}

void Param::set_real(const std::string &key, kids_real val) {
    if (impl_t == JSON) set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
    if (impl_t == TOML) set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val);
}

void Param::set_string(const std::string &key, std::string val) {
    if (impl_t == JSON) set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
    if (impl_t == TOML) set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val);
}

void Param::set_bool_ifndef(const std::string &key, bool val) {
    if (impl_t == JSON) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val); }
    }
    if (impl_t == TOML) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val); }
    }
}

void Param::set_int_ifndef(const std::string &key, int val) {
    if (impl_t == JSON) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val); }
    }
    if (impl_t == TOML) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val); }
    }
}

void Param::set_real_ifndef(const std::string &key, kids_real val) {
    if (impl_t == JSON) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val); }
    }
    if (impl_t == TOML) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val); }
    }
}

void Param::set_string_ifndef(const std::string &key, std::string val) {
    if (impl_t == JSON) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val); }
    }
    if (impl_t == TOML) {
        if (!has_key(key)) { set_internal(std::dynamic_pointer_cast<TomlObject>(obj)->data, key, val); }
    }
}


/**
 * @brief get parameter
 * @param json : DataObject object
 * @param key : key of the parameter
 * @param loc : location tracer (always leaves it `LOC()`)
 * @param qdim : physics dimension for conversion in double in AU unit
 * @param default_value: default value is avialable
 * @return T type
 */
template <class DType, typename T, bool Required = false>
static T get_internal(DType &data, const std::string &key, const std::string &loc,  //
                      const phys::dimension7 &qdim, const T &default_value) {
    auto        ipos       = key.find_first_of(".");
    std::string key1       = (ipos == std::string::npos) ? key : key.substr(0, ipos);
    bool        key1_exist = apply_internal<DType>(data, key1, (bool (*)(const DType &)) internal::is_true);

    if (!key1_exist && Required) {
        throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Required"));
        return T();
    }
    if (!key1_exist && !Required) {
        std::cerr << utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Default=", default_value, "\n");
        return default_value;
    }
    if (ipos != std::string::npos) {
        std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
        if (internal::is_object(data)) {
            return get_internal<DType, T, Required>(internal::child(data, key1), key2, loc, qdim, default_value);
        } else if (internal::is_array(data)) {
            int idx1 = stoi(key1);
            if (idx1 >= 0 && idx1 < data.size())
                return get_internal<DType, T, Required>(data[idx1], key2, loc, qdim, default_value);
        } else {
            throw kids_error("type error");
        }
    }

    // the case find the key
    auto &finaldata = internal::child(data, key);
    if (internal::is_string(finaldata)) {
        if (std::is_same<T, std::string>::value) {
            return internal::as_type<T>(finaldata);
        } else if (std::is_same<T, double>::value) {
            phys::uval        uv   = phys::us::parse(internal::as_string(finaldata));
            double            qval = phys::au::as(qdim, uv);
            T                 q;
            std::stringstream ss;
            ss << std::setiosflags(std::ios::scientific) << std::setprecision(32) << qval;
            ss >> q;
            return q;
        }
    } else if (internal::is_bool(finaldata)) {
        T q;
        if (std::is_same<T, bool>::value) q = internal::as_type<T>(finaldata);
        return q;
    } else if (internal::is_real(finaldata)) {
        T q;
        if (std::is_same<T, double>::value) q = internal::as_real(finaldata);
        return q;
    } else if (internal::is_int(finaldata)) {
        T q;
        if (std::is_same<T, double>::value) q = internal::as_real(finaldata);
        if (std::is_same<T, int>::value) q = internal::as_type<T>(finaldata);
        return q;
    }
    // cannot be adapted to existing conversions
    throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ Data{", internal::dump(finaldata), "}",
                                   " : Converting fatal"));
    return T();
}

template <class DType, typename T, bool Required = false>
T get(DType &data, const std::vector<std::string> &keys, const std::string &loc,  //
      const phys::dimension7 &qdim, const T &default_value) {
    for (auto &key : keys) {
        if (apply_internal<DType>(data, key, (bool (*)(const DType &)) internal::is_true))
            return get_internal<DType, T, Required>(data, key, loc, qdim, default_value);
    }
    if (!Required) return default_value;
    throw kids_error(utils::concat(loc, ": Cannot get parameter for! key = ", keys[0]));
    return T();
}

/// new implementation
/// @{
kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc,
                          const kids_bool &default_value) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_bool, false>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, default_value);
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_bool, false>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, default_value);
    return kids_bool{};
}
kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_bool, true>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, kids_bool());
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_bool, true>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, kids_bool());
    return kids_bool{};
}
kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc, const kids_int &default_value) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_int, false>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, default_value);
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_int, false>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, default_value);
    return kids_int{};
}
kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_int, true>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, kids_int());
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_int, true>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, kids_int());
    return kids_int{};
}
std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc,
                              const std::string &default_value) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, std::string, false>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, default_value);
    if (impl_t == TOML)
        return get<TomlObject::DataType, std::string, false>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, default_value);
    return std::string{};
}
std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, std::string, true>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, std::string());
    if (impl_t == TOML)
        return get<TomlObject::DataType, std::string, true>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, std::string());
    return std::string{};
}

kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc, const phys::dimension7 &qdim,
                          const kids_real &default_value) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_real, false>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, qdim, default_value);
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_real, false>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, qdim, default_value);
    return kids_real{};
}
kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc,
                          const kids_real &default_value) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_real, false>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, default_value);
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_real, false>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, default_value);
    return kids_real{};
}
kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc) {
    if (impl_t == JSON)
        return get<JsonObject::DataType, kids_real, true>(  //
            std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d, kids_real());
    if (impl_t == TOML)
        return get<TomlObject::DataType, kids_real, true>(  //
            std::dynamic_pointer_cast<TomlObject>(obj)->data, keys, loc, phys::none_d, kids_real());
    return kids_real{};
}
/// @}

std::string Param::repr() {
    if (impl_t == JSON) return internal::dump(std::dynamic_pointer_cast<JsonObject>(obj)->data);
    if (impl_t == TOML) return internal::dump(std::dynamic_pointer_cast<TomlObject>(obj)->data);
    return std::string{};
}

// static bool apply_internal(const JsonObject::DataType &data, const std::string &key,
//                            std::function<bool(const JsonObject::DataType &data)> func) {
//     auto        ipos = key.find_first_of(".");
//     std::string key1 = (ipos == std::string::npos) ? key : key.substr(0, ipos);
//     std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
//     if (ipos != std::string::npos) {
//         if (data.is_object()) {
//             if (data.count(key1) == 0) return false;
//             return apply_internal(data[key1], key2, func);
//         } else if (data.is_array()) {
//             int idx1 = stoi(key1);
//             if (idx1 >= 0 && idx1 < data.size()) return apply_internal(data[idx1], key2, func);
//             return false;
//         } else {
//             return false;
//         }
//     }
//     if (data.is_object()) {
//         if (data.count(key1) == 0) return false;
//         return func(data[key]);
//     } else if (data.is_array()) {
//         int idx1 = stoi(key1);
//         if (idx1 >= 0 && idx1 < data.size()) return func(data[idx1]);
//         return false;
//     }
//     return false;
// }

// bool Param::has_key(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return true; });
// }

// bool Param::is_object(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_object(); });
// }

// bool Param::is_array(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_array(); });
// }

// bool Param::is_bool(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_bool(); });
// }

// bool Param::is_int(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_integer(); });
// }

// bool Param::is_real(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_float(); });
// }

// bool Param::is_string(const std::string &key) {
//     return apply_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key,
//                           [](const JsonObject::DataType &data) -> bool { return data.is_string(); });
// }

// template <class T>
// static void set_internal(JsonObject::DataType &data, const std::string &key, const T &val) {
//     auto        ipos = key.find_first_of(".");
//     std::string key1 = (ipos == std::string::npos) ? key : key.substr(0, ipos);
//     std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
//     if (ipos != std::string::npos) {
//         if (data.is_array()) {
//             int idx1 = stoi(key1);
//             if (idx1 >= 0 && idx1 < data.size()) return set_internal(data[idx1], key2, val);
//         }
//         return set_internal(data[key1], key2, val);
//     }
//     data[key] = val;
// }

// void Param::set_bool(const std::string &key, bool val) {
//     set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
// }

// void Param::set_int(const std::string &key, int val) {
//     set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
// }

// void Param::set_real(const std::string &key, kids_real val) {
//     set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
// }

// void Param::set_string(const std::string &key, std::string val) {
//     set_internal(std::dynamic_pointer_cast<JsonObject>(obj)->data, key, val);
// }

// /**
//  * @brief get parameter
//  * @param json : DataObject object
//  * @param key : key of the parameter
//  * @param loc : location tracer (always leaves it `LOC()`)
//  * @param qdim : physics dimension for conversion in double in AU unit
//  * @param default_value: default value is avialable
//  * @return T type
//  */
// template <typename T, bool Required = false>
// static T get_internal(const JsonObject::DataType &data, const std::string &key, const std::string &loc,  //
//                       const phys::dimension7 &qdim, const T &default_value) {
//     auto        ipos       = key.find_first_of(".");
//     std::string key1       = (ipos == std::string::npos) ? key : key.substr(0, ipos);
//     bool        key1_exist =  //
//         apply_internal(data, key1, [](const JsonObject::DataType &data) -> bool { return true; });

//     if (!key1_exist && Required) {
//         throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Required"));
//         return T();
//     }
//     if (!key1_exist && !Required) {
//         std::cerr << utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ : Default=", default_value, "\n");
//         return default_value;
//     }
//     if (ipos != std::string::npos) {
//         std::string key2 = (ipos == std::string::npos) ? "" : key.substr(ipos + 1, key.size());
//         if (data.is_object()) {
//             return get_internal<T, Required>(data[key1], key2, loc, qdim, default_value);
//         } else if (data.is_array()) {
//             int idx1 = stoi(key1);
//             if (idx1 >= 0 && idx1 < data.size())
//                 return get_internal<T, Required>(data[idx1], key2, loc, qdim, default_value);
//         } else {
//             throw kids_error("type error");
//         }
//     }

//     // the case find the key
//     switch (data[key].type()) {
//         case configor::config_value_type::string: {
//             if (std::is_same<T, std::string>::value) {
//                 return data[key].get<T>();
//             } else if (std::is_same<T, double>::value) {
//                 phys::uval        uv   = phys::us::parse(data[key].as_string());
//                 double            qval = phys::au::as(qdim, uv);
//                 T                 q;
//                 std::stringstream ss;
//                 ss << std::setiosflags(std::ios::scientific) << std::setprecision(32) << qval;  // !stupid
//                 ss >> q;
//                 return q;
//             }
//             break;
//         }
//         case configor::config_value_type::boolean: {
//             T q;
//             if (std::is_same<T, bool>::value) q = data[key].get<T>();
//             return q;
//             break;
//         }
//         case configor::config_value_type::number_float: {
//             T q;
//             if (std::is_same<T, double>::value) q = data[key].as_float();
//             return q;
//             break;
//         }
//         case configor::config_value_type::number_integer: {
//             T q;
//             if (std::is_same<T, double>::value) q = data[key].as_float();
//             if (std::is_same<T, int>::value) q = data[key].get<T>();
//             return q;
//             break;
//         }
//     }
//     // cannot be adapted to existing conversions
//     throw kids_error(utils::concat(loc, " Type<", as_str<T>(), "> Key/", key, "/ Data{", data[key].dump(4, ' '), "}",
//                                    " : Converting fatal"));
//     return T();
// }

// template <typename T, bool Required = false>
// T get(const JsonObject::DataType &data, const std::vector<std::string> &keys, const std::string &loc,  //
//       const phys::dimension7 &qdim, const T &default_value) {
//     for (auto &key : keys) {
//         if (apply_internal(data, key, [](const JsonObject::DataType &data) -> bool { return true; }))
//             return get_internal<T, Required>(data, key, loc, qdim, default_value);
//     }
//     if (!Required) return default_value;
//     throw kids_error(utils::concat(loc, ": Cannot get parameter for! key = ", keys[0]));
//     return T();
// }

// /// new implementation
// /// @{
// kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc,
//                           const kids_bool &default_value) {
//     return get<kids_bool, false>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//                                  default_value);
// }
// kids_bool Param::get_bool(const std::vector<std::string> &keys, const std::string &loc) {
//     return get<kids_bool, true>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//     kids_bool());
// }

// kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc, const kids_int &default_value)
// {
//     return get<kids_int, false>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//                                 default_value);
// }
// kids_int Param::get_int(const std::vector<std::string> &keys, const std::string &loc) {
//     return get<kids_int, true>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//     kids_int());
// }

// std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc,
//                               const std::string &default_value) {
//     return get<std::string, false>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//                                    default_value);
// }
// std::string Param::get_string(const std::vector<std::string> &keys, const std::string &loc) {
//     return get<std::string, true>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//                                   std::string());
// }

// kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc, const phys::dimension7 &qdim,
//                           const kids_real &default_value) {
//     return get<kids_real, false>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, qdim, default_value);
// }
// kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc,
//                           const kids_real &default_value) {
//     return get<kids_real, false>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//                                  default_value);
// }
// kids_real Param::get_real(const std::vector<std::string> &keys, const std::string &loc) {
//     return get<kids_real, true>(std::dynamic_pointer_cast<JsonObject>(obj)->data, keys, loc, phys::none_d,
//     kids_real());
// }
// /// @}

// std::string Param::repr() { return std::dynamic_pointer_cast<JsonObject>(obj)->data.dump(4, ' '); }


// #endif  // PARAM_USE_JSON

};  // namespace PROJECT_NS

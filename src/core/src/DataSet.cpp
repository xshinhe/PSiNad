#include "psnd/DataSet.h"

#include <algorithm>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

#include <cstdint>
#include <type_traits>

#include "psnd/Exception.h"
#include "psnd/Shape.h"
#include "psnd/Types.h"
#include "psnd/Variable.h"
#include "psnd/concat.h"

namespace PROJECT_NS {

DataSet::DataSet() {
    _type = psnd_dataset_type;
    _data = std::shared_ptr<DataType>(new DataType());
}

template <typename T>
T* DataSet::def(const std::string& key, Shape S, const std::string& info) {
    DataSetKeyParser          kh    = DataSetKeyParser(key);
    std::shared_ptr<DataType> d_ptr = _data;

    DataSet* currentNode = this;
    for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
        auto& node = (*d_ptr)[kh.terms[i]];
        if (!node) node = std::shared_ptr<DataSet>(new DataSet());

        currentNode = static_cast<DataSet*>(node.get());
        d_ptr       = currentNode->_data;
    }

    std::shared_ptr<Node>& leaf_node = (*d_ptr)[kh.terms.back()];
    if (!leaf_node) {
        // std::cout << "using size for key=" << key << " (" << info << ") : " << S.to_string() << std::endl;
        leaf_node = std::shared_ptr<Tensor<T>>(new Tensor<T>(S, info));
    }

    if (leaf_node->type() != as_enum<T>()) {  //
        std::cout << LOC() << key << "\n";
        throw std::runtime_error("doubly conflicted definition!");
    }
    auto leaf_node_ts = static_cast<Tensor<T>*>(leaf_node.get());
    if (S.size() != leaf_node_ts->size()) {  //
        std::cout << LOC() << key << "\n";
        std::cout << LOC() << S.size() << "\n";
        std::cout << LOC() << leaf_node_ts->size() << "\n";
        std::cout << LOC() << repr() << "\n";
        throw std::runtime_error(utils::concat("doubly conflicted definition!", key));
    }
    return leaf_node_ts->data();
}

template <typename T>
span<T> DataSet::static_def(DataSet& DS, const VARIABLE<T>& var, const span<T>& arr_in, bool allow_diff) {
    span<T> arr(DS.def<T>(var.name(), var.shape(), var.doc()), var.shape().size());
    if (arr_in.size() != 0) {
        if (!allow_diff && arr_in.size() != arr.size()) throw psnd_error("mismatched size when copy dataset data");
        int min_size = std::min(arr.size(), arr_in.size());
        for (int i = 0; i < min_size; ++i) arr[i] = arr_in[i];
    }
    return arr;
}
template <typename T>
span<T> DataSet::static_def(DataSet& DS, const VARIABLE<T>& var, const std::string& ds_file, bool allow_diff) {
    span<T>       arr(DS.def<T>(var.name(), var.shape(), var.doc()), var.shape().size());
    bool          find = false;
    std::string   eachline;
    std::ifstream ifs(ds_file);
    if (!ifs.good()) throw psnd_error(utils::concat("cannot open dataset file from: ", ds_file));
    while (getline(ifs, eachline)) {
        if (eachline == var.name()) {
            getline(ifs, eachline);
            std::string       typeflag;
            std::size_t       vsize;
            std::stringstream ss(eachline);
            ss >> typeflag >> vsize;
            if (!allow_diff && typeflag != as_str<T>()) throw psnd_error("mismatched type");
            if (!allow_diff && vsize != arr.size()) throw psnd_error("mismatched size");
            int min_size = std::min(vsize, arr.size());
            for (int i = 0; i < min_size; ++i) ifs >> arr[i];
            find = true;
        }
    }
    if (!find) throw psnd_error(utils::concat("cannot fetch values in dataset file: ", ds_file));
    return arr;
}
template <typename T>
span<T> DataSet::static_def(DataSet& DS, const VARIABLE<T>& var, std::shared_ptr<DataSet> ds_ptr, bool allow_diff) {
    span<T>    arr(DS.def<T>(var.name(), var.shape(), var.doc()), var.shape().size());
    psnd_dtype itype;
    void*      idata;
    Shape*     ishape;
    std::tie(itype, idata, ishape) = ds_ptr->obtain(var.name());
    if (idata == nullptr) throw psnd_error("ds_ptr has no data");
    if (!allow_diff && itype != as_enum<T>()) throw psnd_error("ds_ptr has mismatch type");
    if (!allow_diff && ishape->size() != var.shape().size()) psnd_error("ds_ptr has mismatch size");
    int min_size = std::min(ishape->size(), var.shape().size());
    for (int i = 0; i < min_size; ++i) arr[i] = ((T*) idata)[i];
    return arr;
}

span<psnd_int> DataSet::def(const VARIABLE<psnd_int>& var, const span<psnd_int>& arr_in) {
    return static_def<psnd_int>(*this, var, arr_in);
}
span<psnd_real> DataSet::def(const VARIABLE<psnd_real>& var, const span<psnd_real>& arr_in) {
    return static_def<psnd_real>(*this, var, arr_in);
}
span<psnd_complex> DataSet::def(const VARIABLE<psnd_complex>& var, const span<psnd_complex>& arr_in) {
    return static_def<psnd_complex>(*this, var, arr_in);
}
span<psnd_int> DataSet::def(const VARIABLE<psnd_int>& var, const std::string& ds_file) {
    return static_def<psnd_int>(*this, var, ds_file);
}
span<psnd_real> DataSet::def(const VARIABLE<psnd_real>& var, const std::string& ds_file) {
    return static_def<psnd_real>(*this, var, ds_file);
}
span<psnd_complex> DataSet::def(const VARIABLE<psnd_complex>& var, const std::string& ds_file) {
    return static_def<psnd_complex>(*this, var, ds_file);
}
span<psnd_int> DataSet::def(const VARIABLE<psnd_int>& var, std::shared_ptr<DataSet> ds_ptr) {
    return static_def<psnd_int>(*this, var, ds_ptr);
}
span<psnd_real> DataSet::def(const VARIABLE<psnd_real>& var, std::shared_ptr<DataSet> ds_ptr) {
    return static_def<psnd_real>(*this, var, ds_ptr);
}
span<psnd_complex> DataSet::def(const VARIABLE<psnd_complex>& var, std::shared_ptr<DataSet> ds_ptr) {
    return static_def<psnd_complex>(*this, var, ds_ptr);
}
void DataSet::def(std::shared_ptr<DataSet> ds_ptr) {
    std::istringstream iss(ds_ptr->repr());
    load(iss);
    return;
};


psnd_int*  DataSet::def_get_pointer(VARIABLE<psnd_int>& var) { return def_int(var.name(), var.shape(), var.doc()); }
psnd_real* DataSet::def_get_pointer(VARIABLE<psnd_real>& var) { return def_real(var.name(), var.shape(), var.doc()); }
psnd_complex* DataSet::def_get_pointer(VARIABLE<psnd_complex>& var) {
    return def_complex(var.name(), var.shape(), var.doc());
}

psnd_int* DataSet::def_int(const std::string& key, Shape S, const std::string& info) {
    return def<psnd_int>(key, S, info);
}
psnd_int* DataSet::def_int(const std::string& key, psnd_int* arr_in, Shape S, const std::string& info) {
    psnd_int* arr = def_int(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
psnd_int* DataSet::def_int(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == psnd_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<psnd_int>*>(inode);
    return def_int(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_int(const std::string& key, Shape S, const std::string& info) {
    def_int(key, S, info);
    return *this;
}

psnd_real* DataSet::def_real(const std::string& key, Shape S, const std::string& info) {
    return def<psnd_real>(key, S, info);
}
psnd_real* DataSet::def_real_replace(const std::string& key, Shape S, const std::string& info) {
    if (!haskey(key)) return def<psnd_real>(key, S, info);
    // return def<psnd_real>(key, S, info);

    std::cout << LOC() << S.to_string() << "\n";

    auto       old_node = obtain(key);
    int        size_min = std::min(S.size(), std::get<2>(old_node)->size());
    psnd_real* old_arr  = (psnd_real*) std::get<1>(old_node);

    std::string tmp_key = utils::concat("tmpr.", key);
    psnd_real*  arr     = def_real(tmp_key, S, info);

    std::cout << LOC() << key << "\n";
    // std::cout << LOC() << tmp_key << "\n";
    // std::cout << LOC() << "min = " << size_min << "\n";
    for (int i = 0; i < size_min; ++i) arr[i] = old_arr[i];
    _undef(key);
    auto res = def_real(key, tmp_key, info);
    _undef(tmp_key);
    // std::cout << LOC() << key << "\n";
    return res;
}
psnd_real* DataSet::def_real(const std::string& key, psnd_real* arr_in, Shape S, const std::string& info) {
    psnd_real* arr = def_real(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
psnd_real* DataSet::def_real(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == psnd_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<psnd_real>*>(inode);
    return def_real(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_real(const std::string& key, Shape S, const std::string& info) {
    def_real(key, S, info);
    return *this;
}

psnd_complex* DataSet::def_complex(const std::string& key, Shape S, const std::string& info) {
    return def<psnd_complex>(key, S, info);
}
psnd_complex* DataSet::def_complex_replace(const std::string& key, Shape S, const std::string& info) {
    if (!haskey(key)) return def<psnd_complex>(key, S, info);
    // return def<psnd_complex>(key, S, info);

    auto          old_node = obtain(key);
    int           size_min = std::min(S.size(), std::get<2>(old_node)->size());
    psnd_complex* old_arr  = (psnd_complex*) std::get<1>(old_node);

    std::string   tmp_key = utils::concat("tmpc.", key);
    psnd_complex* arr     = def_complex(tmp_key, S, info);

    std::cout << LOC() << key << "\n";
    // std::cout << LOC() << tmp_key << "\n";
    // std::cout << LOC() << "min = " << size_min << "\n";
    for (int i = 0; i < size_min; ++i) arr[i] = old_arr[i];
    _undef(key);
    auto res = def_complex(key, tmp_key, info);
    _undef(tmp_key);
    // std::cout << LOC() << key << "\n";

    return res;
}
psnd_complex* DataSet::def_complex(const std::string& key, psnd_complex* arr_in, Shape S, const std::string& info) {
    psnd_complex* arr = def_complex(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
psnd_complex* DataSet::def_complex(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == psnd_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<psnd_complex>*>(inode);
    return def_complex(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_complex(const std::string& key, Shape S, const std::string& info) {
    def_complex(key, S, info);
    return *this;
}

DataSet& DataSet::_def(const std::string& key, const std::string& key_in, const std::string& info) {
    auto leaf_node = node(key_in);
    switch (leaf_node->type()) {
        case psnd_int_type:
            def_int(key, key_in, info);
            break;
        case psnd_real_type:
            def_real(key, key_in, info);
            break;
        case psnd_complex_type:
            def_complex(key, key_in, info);
            break;
        default:
            break;
    }
    return *this;
}

DataSet& DataSet::_undef(const std::string& key) {
    DataSetKeyParser          kh    = DataSetKeyParser(key);
    std::shared_ptr<DataType> d_ptr = _data;

    DataSet* currentNode = this;
    for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
        auto& node = (*d_ptr)[kh.terms[i]];
        if (!node) return *this;

        currentNode = static_cast<DataSet*>(node.get());
        d_ptr       = currentNode->_data;
    }
    auto it = d_ptr->find(kh.terms.back());
    if (it != d_ptr->end()) {
        std::cout << LOC() << key << " help!!!\n";
        // it->second.reset();
        // d_ptr->erase(it);  // @TODO
    }
    return *this;
}

std::tuple<psnd_dtype, void*, Shape*> DataSet::obtain(const std::string& key) {
    auto&& leaf_node = node(key);
    switch (leaf_node->type()) {
        case psnd_int_type: {
            auto&& conv_node = static_cast<Tensor<psnd_int>*>(leaf_node);
            return std::make_tuple(psnd_int_type, conv_node->data(), &(conv_node->shape()));
            break;
        }
        case psnd_real_type: {
            auto&& conv_node = static_cast<Tensor<psnd_real>*>(leaf_node);
            return std::make_tuple(psnd_real_type, conv_node->data(), &(conv_node->shape()));
            break;
        }
        case psnd_complex_type: {
            auto&& conv_node = static_cast<Tensor<psnd_complex>*>(leaf_node);
            return std::make_tuple(psnd_complex_type, conv_node->data(), &(conv_node->shape()));
            break;
        }
        default: {
            throw std::runtime_error("bad obtain!");
            break;
        }
    }
}

bool DataSet::haskey(const std::string& key) {
    DataSetKeyParser          kh    = DataSetKeyParser(key);
    std::shared_ptr<DataType> d_ptr = _data;

    DataSet* currentNode = this;
    for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
        auto& node = (*d_ptr)[kh.terms[i]];
        if (!node) return false;
        currentNode = static_cast<DataSet*>(node.get());
        d_ptr       = currentNode->_data;
    }

    auto& leaf_node = (*d_ptr)[kh.terms.back()];
    if (!leaf_node) return false;
    return true;
}

Node* DataSet::node(const std::string& key) {
    DataSetKeyParser          kh    = DataSetKeyParser(key);
    std::shared_ptr<DataType> d_ptr = _data;

    DataSet* currentNode = this;
    for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
        auto& node = (*d_ptr)[kh.terms[i]];
        if (!node) {
            std::cout << "[ERROR] node() failed at intermediate term: '" << kh.terms[i] << "' for key: '" << key << "'"
                      << std::endl;
            std::cout << "[ERROR] This suggests the key is being accessed (not defined) before it exists" << std::endl;
            throw std::runtime_error(std::string{key} + " : access undefined key!");
        }

        currentNode = static_cast<DataSet*>(node.get());
        d_ptr       = currentNode->_data;
    }

    auto& leaf_node = (*d_ptr)[kh.terms.back()];
    if (!leaf_node) throw std::runtime_error(std::string{key} + " : access undefined key!");
    return leaf_node.get();
}

DataSet* DataSet::at(const std::string& key) {
    auto leaf_node = node(key);
    if (leaf_node->type() == psnd_dataset_type) {
        return static_cast<DataSet*>(leaf_node);
    } else {
        throw std::runtime_error("bad conversion!");
    }
    return nullptr;
}

std::string DataSet::help(const std::string& name) {
    std::ostringstream                          os;
    std::vector<std::tuple<std::string, Node*>> stack;

    Node* inode = this;
    if (name != "") inode = node(name);

    stack.push_back(std::make_tuple("", inode));
    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back();
        stack.pop_back();

        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;

        for (auto& i : (*d_ptr)) {
            std::string key   = (parent == "") ? i.first : parent + "." + i.first;
            Node*       inode = i.second.get();
            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                os << key << ":\n\t" << inode->help("") << "\n";
            }
        }
    }
    return os.str();
}

std::string DataSet::repr() {
    std::ostringstream                          os;
    std::vector<std::tuple<std::string, Node*>> stack;

    stack.push_back(std::make_tuple("", this));

    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back();
        stack.pop_back();
        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
        for (auto& i : (*d_ptr)) {
            std::string key = (parent == "") ? i.first : parent + "." + i.first;
            if (!i.second) continue;
            Node* inode = i.second.get();

            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                os << key << "\n" << inode->repr() << "\n\n";
            }
        }
    }
    return os.str();
}

void DataSet::dump_match(std::ostream& os, const std::string& prefix) {
    std::vector<std::tuple<std::string, Node*>> stack;
    stack.push_back(std::make_tuple("", this));

    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back();
        stack.pop_back();
        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
        for (auto& i : (*d_ptr)) {
            std::string key = (parent == "") ? i.first : parent + "." + i.first;
            if (!i.second) continue;
            Node* inode = i.second.get();

            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                auto ipos = key.find(prefix);
                if (ipos != std::string::npos && ipos == 0) {  //
                    os << key << "\n" << inode->repr() << "\n\n";
                }
            }
        }
    }
}

void DataSet::dump(std::ostream& os) { os << repr(); }

void DataSet::load(std::istream& is) {
    std::string key, typeflag, eachline;
    std::size_t size;
    while (is >> key) {
        getline(is, eachline);
        getline(is, eachline);
        std::stringstream ss(eachline);
        ss >> typeflag >> size;
        std::size_t              idim;
        std::vector<std::size_t> dims;
        while (ss >> idim) dims.push_back(idim);
        Shape shtmp(dims);
        if (shtmp.size() != size) throw psnd_error("load ds error");
        if (typeflag.find(as_str<psnd_int>()) != std::string::npos) {
            // nsamp should be carefully checked with Param!!! @bug
            int* ptr = def<int>(key, shtmp);
            for (int i = 0; i < size; ++i) {
                std::string tok;
                is >> tok;
                try {
                    ptr[i] = std::stoi(tok);
                } catch (const std::invalid_argument& e) {
                    throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                }
            }
            // check isnan and isinf
            bool has_nan = std::any_of(ptr, ptr + size, [](psnd_int v) { return false; });
            bool has_inf = std::any_of(ptr, ptr + size, [](psnd_int v) { return false; });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        } else if (typeflag.find(as_str<psnd_real>()) != std::string::npos) {
            psnd_real* ptr = def<psnd_real>(key, shtmp);
            // std::cout << LOC() << is.rdbuf() << std::endl;
            for (int i = 0; i < size; ++i) {
                // is >> ptr[i]; 
                std::string tok;
                is >> tok;
                try {
                    ptr[i] = std::stod(tok);
                } catch (const std::invalid_argument& e) {
                    throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                } 
            }
            // check isnan
            bool has_nan = std::any_of(ptr, ptr + size, [](psnd_real v) { return std::isnan(v); });
            bool has_inf = std::any_of(ptr, ptr + size, [](psnd_real v) { return std::isinf(v); });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        } else if (typeflag.find(as_str<psnd_complex>()) != std::string::npos) {
            psnd_complex* ptr = def<psnd_complex>(key, shtmp);
            for (int i = 0; i < size; ++i) {
                std::string tok;
                is >> tok;
                // string to complex
                // 假设 tok 为 "(1.472428460107e-01,9.891003711953e-01)"
                size_t pos1 = tok.find('(');
                size_t pos2 = tok.find(',');
                size_t pos3 = tok.find(')');
                try{
                    psnd_real real = std::stod(tok.substr(pos1 + 1, pos2 - pos1 - 1));
                    psnd_real imag = std::stod(tok.substr(pos2 + 1, pos3 - pos2 - 1));
                    ptr[i] = std::complex<psnd_real>{real, imag};
                } catch (const std::invalid_argument& e) {
                    throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                } 
            }
            // check isnan and inf
            bool has_nan = std::any_of(ptr, ptr + size, [](psnd_complex v) {
                return std::isnan(v.real()) || std::isnan(v.imag());
            });
            bool has_inf = std::any_of(ptr, ptr + size, [](psnd_complex v) {
                return std::isinf(v.real()) || std::isinf(v.imag());
            });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        }
    }
}

void DataSet::load_reframe(std::istream& is, std::size_t nsamp) {
    std::string key, typeflag, eachline;
    int         size;
    
    while (is >> key) {
        std::cout << LOC() << "loading key = " << key << "\n";
        bool multiframe = (key[0] == '_');
        getline(is, eachline);
        getline(is, eachline);
        std::cout << LOC() << "header line = " << eachline << "\n";
        std::stringstream ss(eachline);
        ss >> typeflag >> size;
        std::size_t              idim;
        std::vector<std::size_t> dims;
        while (ss >> idim) dims.push_back(idim);
        if (multiframe) dims[0] = nsamp;  // update shape with nsamp
        Shape shtmp(dims);

        std::cout << LOC() << "expected size = " << shtmp.size() << " actual size = " << size << "\n";

        int min_size = std::min(size, shtmp.size());
        if (typeflag.find(as_str<psnd_int>()) != std::string::npos) {
            int* ptr = def<int>(key, shtmp);
            int tmpi;
            if (key == "control.nsamp") {
                for (int i = 0; i < size; ++i) is >> tmpi;
            } else {
                for (int i = 0; i < min_size; ++i) {
                    std::string tok;
                    if (!(is >> tok)) {
                        throw psnd_error(utils::concat("load_reframe: unexpected EOF while reading values for ", key));
                    }
                    try {
                        ptr[i] = std::stoi(tok);
                    } catch (const std::invalid_argument& e) {
                        throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                    }
                }
            }
            // check isnan and isinf (for int, usually not needed, but keep for consistency)
            bool has_nan = std::any_of(ptr, ptr + min_size, [](psnd_int v) { return false; });
            bool has_inf = std::any_of(ptr, ptr + min_size, [](psnd_int v) { return false; });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        } else if (typeflag.find(as_str<psnd_real>()) != std::string::npos) {
            psnd_real* ptr = def<psnd_real>(key, shtmp);
            for (int i = 0; i < min_size; ++i) {
                std::string tok;
                if (!(is >> tok)) {
                    throw psnd_error(utils::concat("load_reframe: unexpected EOF while reading values for ", key));
                }
                try {
                    ptr[i] = std::stod(tok);
                } catch (const std::invalid_argument& e) {
                    throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                }
            }
            bool has_nan = std::any_of(ptr, ptr + min_size, [](psnd_real v) { return std::isnan(v); });
            bool has_inf = std::any_of(ptr, ptr + min_size, [](psnd_real v) { return std::isinf(v); });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        } else if (typeflag.find(as_str<psnd_complex>()) != std::string::npos) {
            psnd_complex* ptr = def<psnd_complex>(key, shtmp);
            for (int i = 0; i < min_size; ++i) {
                std::string tok;
                if (!(is >> tok)) {
                    throw psnd_error(utils::concat("load_reframe: unexpected EOF while reading values for ", key));
                }
                // string to complex, e.g. "(1.472428460107e-01,9.891003711953e-01)"
                size_t pos1 = tok.find('(');
                size_t pos2 = tok.find(',');
                size_t pos3 = tok.find(')');
                try {
                    if (pos1 == std::string::npos || pos2 == std::string::npos || pos3 == std::string::npos) {
                        throw std::invalid_argument("bad complex format");
                    }
                    psnd_real real = std::stod(tok.substr(pos1 + 1, pos2 - pos1 - 1));
                    psnd_real imag = std::stod(tok.substr(pos2 + 1, pos3 - pos2 - 1));
                    if (std::isinf(real) || std::isinf(imag)) {
                        throw std::runtime_error("Infinite value encountered in input for key: " + key);
                    }
                    ptr[i] = std::complex<psnd_real>{real, imag};
                } catch (const std::invalid_argument& e) {
                    throw psnd_error("Invalid input for key: " + key + ", value: " + tok + ", index: " + std::to_string(i));
                }
            }
            bool has_nan = std::any_of(ptr, ptr + min_size, [](psnd_complex v) {
                return std::isnan(v.real()) || std::isnan(v.imag());
            });
            bool has_inf = std::any_of(ptr, ptr + min_size, [](psnd_complex v) {
                return std::isinf(v.real()) || std::isinf(v.imag());
            });
            if (has_nan || has_inf) {
                throw psnd_error("NaN or Inf value encountered in input for key: " + key);
            }
        }
    }
    std::cout << LOC() << "load_reframe done\n";
}


/*
----------------------------------------------NEW FEATURE----------------------------------------------------
have not been totally reviewed and tested yet,

using binary for saving and loading datasets could be implemented here for efficiency. hclu 251027
*/

// Estimate binary size (bytes) if dataset serialized with the format below.
std::size_t DataSet::estimate_binary_size() const {
    std::size_t total = 0;
    // 4 bytes header: number of entries (uint32_t)
    total += sizeof(uint32_t);
    std::vector<std::tuple<std::string, Node*>> stack;
    stack.push_back(std::make_tuple(std::string(""), const_cast<DataSet*>(this)));
    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back();
        stack.pop_back();
        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
        for (auto &i : (*d_ptr)) {
            if (!i.second) continue;
            std::string key = (parent == "") ? i.first : parent + "." + i.first;
            Node* inode = i.second.get();
            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                // key length(4) + key bytes + type(1) + rank(4) + dims(uint64_t * rank) + data bytes
                total += sizeof(uint32_t) + key.size() + sizeof(uint8_t) + sizeof(uint32_t);
                // shape
                Shape *sh = nullptr;
                switch (inode->type()) {
                    case psnd_int_type: sh = &(static_cast<Tensor<psnd_int>*>(inode)->shape()); break;
                    case psnd_real_type: sh = &(static_cast<Tensor<psnd_real>*>(inode)->shape()); break;
                    case psnd_complex_type: sh = &(static_cast<Tensor<psnd_complex>*>(inode)->shape()); break;
                    default: break;
                }
                if (sh) {
                    total += sizeof(uint64_t) * sh->rank();
                    std::size_t cnt = sh->size();
                    switch (inode->type()) {
                        case psnd_int_type: total += cnt * sizeof(psnd_int); break;
                        case psnd_real_type: total += cnt * sizeof(psnd_real); break;
                        case psnd_complex_type: total += cnt * sizeof(psnd_complex); break;
                        default: break;
                    }
                }
            }
        }
    }
    return total;
}

// Save dataset in compact binary format:
// [uint32_t] n_entries
// for each entry:
// [uint32_t] key_len
// [key_len bytes] key (no null)
// [uint8_t] dtype (0=int,1=real,2=complex)
// [uint32_t] rank
// [uint64_t]*rank dims
// [data bytes] raw elements (row-major)
void DataSet::save_binary(const std::string &path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs) throw psnd_error(utils::concat("cannot open binary file for write: ", path));
    // first pass: count entries
    uint32_t count = 0;
    {
        std::vector<std::tuple<std::string, Node*>> stack;
        stack.push_back(std::make_tuple(std::string(""), const_cast<DataSet*>(this)));
        while (!stack.empty()) {
            auto [parent, currentNode] = stack.back(); stack.pop_back();
            std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
            for (auto &i : (*d_ptr)) {
                if (!i.second) continue;
                std::string key = (parent == "") ? i.first : parent + "." + i.first;
                Node* inode = i.second.get();
                if (inode->type() == psnd_dataset_type) {
                    stack.push_back(std::make_tuple(key, inode));
                } else {
                    ++count;
                }
            }
        }
    }
    ofs.write(reinterpret_cast<const char*>(&count), sizeof(count));
    // second pass: write entries
    std::vector<std::tuple<std::string, Node*>> stack;
    stack.push_back(std::make_tuple(std::string(""), const_cast<DataSet*>(this)));
    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back(); stack.pop_back();
        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
        for (auto &i : (*d_ptr)) {
            if (!i.second) continue;
            std::string key = (parent == "") ? i.first : parent + "." + i.first;
            Node* inode = i.second.get();
            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                // key
                uint32_t klen = static_cast<uint32_t>(key.size());
                ofs.write(reinterpret_cast<const char*>(&klen), sizeof(klen));
                ofs.write(key.data(), klen);
                // dtype
                uint8_t dtype = 0;
                if (inode->type() == psnd_real_type) dtype = 1;
                else if (inode->type() == psnd_complex_type) dtype = 2;
                ofs.write(reinterpret_cast<const char*>(&dtype), sizeof(dtype));
                // shape
                Shape *sh = nullptr;
                if (inode->type() == psnd_int_type) sh = &(static_cast<Tensor<psnd_int>*>(inode)->shape());
                else if (inode->type() == psnd_real_type) sh = &(static_cast<Tensor<psnd_real>*>(inode)->shape());
                else if (inode->type() == psnd_complex_type) sh = &(static_cast<Tensor<psnd_complex>*>(inode)->shape());
                uint32_t rank = static_cast<uint32_t>(sh->rank());
                ofs.write(reinterpret_cast<const char*>(&rank), sizeof(rank));
                for (auto d : sh->dims()) {
                    uint64_t ud = static_cast<uint64_t>(d);
                    ofs.write(reinterpret_cast<const char*>(&ud), sizeof(ud));
                }
                // data
                std::size_t cnt = sh->size();
                if (inode->type() == psnd_int_type) {
                    auto ptr = static_cast<Tensor<psnd_int>*>(inode)->data();
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_int));
                } else if (inode->type() == psnd_real_type) {
                    auto ptr = static_cast<Tensor<psnd_real>*>(inode)->data();
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_real));
                } else if (inode->type() == psnd_complex_type) {
                    auto ptr = static_cast<Tensor<psnd_complex>*>(inode)->data();
                    static_assert(std::is_trivially_copyable<psnd_complex>::value, "complex must be trivially copyable");
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_complex));
                }
            }
        }
    }
    ofs.close();
}

void DataSet::save_binary_filter(const std::string &path) const {
    std::ofstream ofs(path, std::ios::binary);
    if (!ofs)
        throw psnd_error(utils::concat("cannot open binary file for write: ", path));

    // ---------- First pass: count entries ----------
    uint32_t count = 0;
    {
        std::vector<std::tuple<std::string, Node*>> stack;
        stack.push_back(std::make_tuple(std::string(""), const_cast<DataSet*>(this)));
        while (!stack.empty()) {
            auto [parent, currentNode] = stack.back(); stack.pop_back();
            std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
            for (auto &i : (*d_ptr)) {
                if (!i.second) continue;
                std::string key = (parent.empty()) ? i.first : parent + "." + i.first;
                Node* inode = i.second.get();
                if (inode->type() == psnd_dataset_type) {
                    stack.push_back(std::make_tuple(key, inode));
                } else {
                    // 🔹 忽略以 "_.1" 或 "_.2" 开头的 key
                    if (key.rfind("_.1", 0) == 0 || key.rfind("_.2", 0) == 0)
                        continue;
                    ++count;
                }
            }
        }
    }

    ofs.write(reinterpret_cast<const char*>(&count), sizeof(count));

    // ---------- Second pass: write entries ----------
    std::vector<std::tuple<std::string, Node*>> stack;
    stack.push_back(std::make_tuple(std::string(""), const_cast<DataSet*>(this)));
    while (!stack.empty()) {
        auto [parent, currentNode] = stack.back(); stack.pop_back();
        std::shared_ptr<DataType> d_ptr = static_cast<DataSet*>(currentNode)->_data;
        for (auto &i : (*d_ptr)) {
            if (!i.second) continue;
            std::string key = (parent.empty()) ? i.first : parent + "." + i.first;
            Node* inode = i.second.get();
            if (inode->type() == psnd_dataset_type) {
                stack.push_back(std::make_tuple(key, inode));
            } else {
                // 🔹 忽略以 "_.1" 或 "_.2" 开头的 key
                if (key.rfind("_.1", 0) == 0 || key.rfind("_.2", 0) == 0)
                    continue;

                // ---------- 写 key ----------
                uint32_t klen = static_cast<uint32_t>(key.size());
                ofs.write(reinterpret_cast<const char*>(&klen), sizeof(klen));
                ofs.write(key.data(), klen);

                // ---------- 写 dtype ----------
                uint8_t dtype = 0;
                if (inode->type() == psnd_real_type)
                    dtype = 1;
                else if (inode->type() == psnd_complex_type)
                    dtype = 2;
                ofs.write(reinterpret_cast<const char*>(&dtype), sizeof(dtype));

                // ---------- 写 shape ----------
                Shape* sh = nullptr;
                if (inode->type() == psnd_int_type)
                    sh = &(static_cast<Tensor<psnd_int>*>(inode)->shape());
                else if (inode->type() == psnd_real_type)
                    sh = &(static_cast<Tensor<psnd_real>*>(inode)->shape());
                else if (inode->type() == psnd_complex_type)
                    sh = &(static_cast<Tensor<psnd_complex>*>(inode)->shape());

                uint32_t rank = static_cast<uint32_t>(sh->rank());
                ofs.write(reinterpret_cast<const char*>(&rank), sizeof(rank));
                for (auto d : sh->dims()) {
                    uint64_t ud = static_cast<uint64_t>(d);
                    ofs.write(reinterpret_cast<const char*>(&ud), sizeof(ud));
                }

                // ---------- 写数据 ----------
                std::size_t cnt = sh->size();
                if (inode->type() == psnd_int_type) {
                    auto ptr = static_cast<Tensor<psnd_int>*>(inode)->data();
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_int));
                } else if (inode->type() == psnd_real_type) {
                    auto ptr = static_cast<Tensor<psnd_real>*>(inode)->data();
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_real));
                } else if (inode->type() == psnd_complex_type) {
                    auto ptr = static_cast<Tensor<psnd_complex>*>(inode)->data();
                    static_assert(std::is_trivially_copyable<psnd_complex>::value, "complex must be trivially copyable");
                    ofs.write(reinterpret_cast<const char*>(ptr), cnt * sizeof(psnd_complex));
                }
            }
        }
    }

    ofs.close();
}



// Load dataset from the binary format above. Existing keys will be overwritten/defined via def<>().
void DataSet::load_binary(const std::string &path) {
    std::ifstream ifs(path, std::ios::binary);
    if (!ifs) throw psnd_error(utils::concat("cannot open binary file for read: ", path));
    uint32_t count = 0;
    ifs.read(reinterpret_cast<char*>(&count), sizeof(count));
    for (uint32_t e = 0; e < count; ++e) {
        uint32_t klen = 0;
        ifs.read(reinterpret_cast<char*>(&klen), sizeof(klen));
        std::string key(klen, '\0');
        ifs.read(&key[0], klen);
        uint8_t dtype = 0;
        ifs.read(reinterpret_cast<char*>(&dtype), sizeof(dtype));
        uint32_t rank = 0;
        ifs.read(reinterpret_cast<char*>(&rank), sizeof(rank));
        std::vector<std::size_t> dims(rank);
        for (uint32_t r = 0; r < rank; ++r) {
            uint64_t ud = 0;
            ifs.read(reinterpret_cast<char*>(&ud), sizeof(ud));
            dims[r] = static_cast<std::size_t>(ud);
        }
        Shape shtmp(dims);
        std::size_t cnt = shtmp.size();
        if (dtype == 0) {
            psnd_int* ptr = def<psnd_int>(key, shtmp);
            ifs.read(reinterpret_cast<char*>(ptr), cnt * sizeof(psnd_int));
        } else if (dtype == 1) {
            psnd_real* ptr = def<psnd_real>(key, shtmp);
            ifs.read(reinterpret_cast<char*>(ptr), cnt * sizeof(psnd_real));
        } else if (dtype == 2) {
            psnd_complex* ptr = def<psnd_complex>(key, shtmp);
            ifs.read(reinterpret_cast<char*>(ptr), cnt * sizeof(psnd_complex));
        } else {
            throw psnd_error("unknown dtype in binary file");
        }
        if (!ifs) throw psnd_error(utils::concat("unexpected EOF or IO error while reading ", path));
    }
    ifs.close();
}


};  // namespace PROJECT_NS

/**
int main() {
    using namespace PROJECT_NS;

    DataSet DS;
    DS.def<int>("0.1", 4);
    DS.def<int>("a.b", 10);
    DS.def<double>("a.c.1", 8);
    DS.def<double>("a.c.d", 8);
    DS.def<double>("a.c.f", Shape({1, 2, 3}));

    DS.dump(std::cout);
    DS._undef("a.c.d");
    DS.def<int>("a.c.d", 3);
    DS.dump(std::cout);

    auto&& DS2 = DS.at("a");
    DS2->dump(std::cout);

    return 0;
}
*/

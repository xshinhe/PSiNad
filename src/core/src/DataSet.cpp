#include "kids/DataSet.h"

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

#include "kids/Exception.h"
#include "kids/Shape.h"
#include "kids/Types.h"
#include "kids/Variable.h"
#include "kids/concat.h"

namespace PROJECT_NS {

DataSet::DataSet() {
    _type = kids_dataset_type;
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
    if (!leaf_node) leaf_node = std::shared_ptr<Tensor<T>>(new Tensor<T>(S, info));

    if (leaf_node->type() != as_enum<T>()) {  //
        throw std::runtime_error("doubly conflicted definition!");
    }
    auto leaf_node_ts = static_cast<Tensor<T>*>(leaf_node.get());
    if (S.size() != leaf_node_ts->size()) {  //
        throw std::runtime_error("doubly conflicted definition!");
    }
    return leaf_node_ts->data();
}

template <typename T>
span<T> DataSet::static_def(DataSet& DS, const VARIABLE<T>& var, const span<T>& arr_in) {
    span<T> arr(DS.def<T>(var.name(), var.shape(), var.doc()), var.shape().size());
    if (arr_in.size() != 0) {
        if (arr_in.size() != arr.size()) throw kids_error("mismatched size when copy dataset data");
        for (int i = 0; i < arr.size(); ++i) arr[i] = arr_in[i];
    }
    return arr;
}

template <typename T>
span<T> DataSet::static_def(DataSet& DS, const VARIABLE<T>& var, const std::string& ds_file) {
    span<T> arr(DS.def<T>(var.name(), var.shape(), var.doc()), var.shape().size());

    bool          find = false;
    std::string   eachline;
    std::ifstream ifs(ds_file);
    if (!ifs.good()) throw kids_error(utils::concat("cannot open dataset file from: ", ds_file));
    while (getline(ifs, eachline)) {
        if (eachline == var.name()) {
            getline(ifs, eachline);
            std::string       typeflag;
            int               vsize;
            std::stringstream ss(eachline);
            ss >> typeflag >> vsize;
            if (typeflag != as_str<T>()) throw kids_error("mismatched type");
            if (vsize != arr.size()) throw kids_error("mismatched size");
            for (int i = 0; i < arr.size(); ++i) ifs >> arr[i];
            find = true;
        }
    }
    if (!find) throw kids_error(utils::concat("cannot fetch values in dataset file: ", ds_file));
    return arr;
}

span<kids_int> DataSet::def(const VARIABLE<kids_int>& var, const span<kids_int>& arr_in) {
    return static_def<kids_int>(*this, var, arr_in);
}
span<kids_real> DataSet::def(const VARIABLE<kids_real>& var, const span<kids_real>& arr_in) {
    return static_def<kids_real>(*this, var, arr_in);
}
span<kids_complex> DataSet::def(const VARIABLE<kids_complex>& var, const span<kids_complex>& arr_in) {
    return static_def<kids_complex>(*this, var, arr_in);
}
span<kids_int> DataSet::def(const VARIABLE<kids_int>& var, const std::string& ds_file) {
    return static_def<kids_int>(*this, var, ds_file);
}
span<kids_real> DataSet::def(const VARIABLE<kids_real>& var, const std::string& ds_file) {
    return static_def<kids_real>(*this, var, ds_file);
}
span<kids_complex> DataSet::def(const VARIABLE<kids_complex>& var, const std::string& ds_file) {
    return static_def<kids_complex>(*this, var, ds_file);
}

kids_int*  DataSet::def_get_pointer(VARIABLE<kids_int>& var) { return def_int(var.name(), var.shape(), var.doc()); }
kids_real* DataSet::def_get_pointer(VARIABLE<kids_real>& var) { return def_real(var.name(), var.shape(), var.doc()); }
kids_complex* DataSet::def_get_pointer(VARIABLE<kids_complex>& var) {
    return def_complex(var.name(), var.shape(), var.doc());
}


kids_int* DataSet::def_int(const std::string& key, Shape S, const std::string& info) {
    return def<kids_int>(key, S, info);
}
kids_int* DataSet::def_int(const std::string& key, kids_int* arr_in, Shape S, const std::string& info) {
    kids_int* arr = def_int(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
kids_int* DataSet::def_int(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == kids_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<kids_int>*>(inode);
    return def_int(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_int(const std::string& key, Shape S, const std::string& info) {
    def_int(key, S, info);
    return *this;
}

kids_real* DataSet::def_real(const std::string& key, Shape S, const std::string& info) {
    return def<kids_real>(key, S, info);
}
kids_real* DataSet::def_real(const std::string& key, kids_real* arr_in, Shape S, const std::string& info) {
    kids_real* arr = def_real(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
kids_real* DataSet::def_real(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == kids_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<kids_real>*>(inode);
    return def_real(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_real(const std::string& key, Shape S, const std::string& info) {
    def_real(key, S, info);
    return *this;
}

kids_complex* DataSet::def_complex(const std::string& key, Shape S, const std::string& info) {
    return def<kids_complex>(key, S, info);
}
kids_complex* DataSet::def_complex(const std::string& key, kids_complex* arr_in, Shape S, const std::string& info) {
    kids_complex* arr = def_complex(key, S, info);
    for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
    return arr;
}
kids_complex* DataSet::def_complex(const std::string& key, const std::string& key_in, const std::string& info) {
    auto inode = node(key_in);
    if (inode->type() == kids_dataset_type) {  //
        throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
    }
    auto inode_ts = static_cast<Tensor<kids_complex>*>(inode);
    return def_complex(key, inode_ts->data(), inode_ts->size(), info);
}
DataSet& DataSet::_def_complex(const std::string& key, Shape S, const std::string& info) {
    def_complex(key, S, info);
    return *this;
}


DataSet& DataSet::_def(const std::string& key, const std::string& key_in, const std::string& info) {
    auto leaf_node = node(key_in);
    switch (leaf_node->type()) {
        case kids_int_type:
            def_int(key, key_in, info);
            break;
        case kids_real_type:
            def_real(key, key_in, info);
            break;
        case kids_complex_type:
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
        std::cout << "\n";
        it->second.reset();
    }
    return *this;
}

std::tuple<kids_dtype, void*, Shape*> DataSet::obtain(const std::string& key) {
    auto&& leaf_node = node(key);
    switch (leaf_node->type()) {
        case kids_int_type: {
            auto&& conv_node = static_cast<Tensor<kids_int>*>(leaf_node);
            return std::make_tuple(kids_int_type, conv_node->data(), &(conv_node->shape()));
            break;
        }
        case kids_real_type: {
            auto&& conv_node = static_cast<Tensor<kids_real>*>(leaf_node);
            return std::make_tuple(kids_real_type, conv_node->data(), &(conv_node->shape()));
            break;
        }
        case kids_complex_type: {
            auto&& conv_node = static_cast<Tensor<kids_complex>*>(leaf_node);
            return std::make_tuple(kids_complex_type, conv_node->data(), &(conv_node->shape()));
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
        if (!node) throw std::runtime_error(std::string{key} + " : access undefined key!");

        currentNode = static_cast<DataSet*>(node.get());
        d_ptr       = currentNode->_data;
    }

    auto& leaf_node = (*d_ptr)[kh.terms.back()];
    if (!leaf_node) throw std::runtime_error(std::string{key} + " : access undefined key!");
    return leaf_node.get();
}

DataSet* DataSet::at(const std::string& key) {
    auto leaf_node = node(key);
    if (leaf_node->type() == kids_dataset_type) {
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
            if (inode->type() == kids_dataset_type) {
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

            if (inode->type() == kids_dataset_type) {
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

            if (inode->type() == kids_dataset_type) {
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
    std::string key, typeflag;
    int         size;
    while (is >> key >> typeflag >> size) {
        if (typeflag == as_str<int>()) {
            int* ptr = def<int>(key, size);
            for (int i = 0; i < size; ++i) is >> ptr[i];
        } else if (typeflag == as_str<kids_real>()) {
            kids_real* ptr = def<kids_real>(key, size);
            for (int i = 0; i < size; ++i) is >> ptr[i];
        } else if (typeflag == as_str<kids_complex>()) {
            kids_complex* ptr = def<kids_complex>(key, size);
            for (int i = 0; i < size; ++i) is >> ptr[i];
        }
    }
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

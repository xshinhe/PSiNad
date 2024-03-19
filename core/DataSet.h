/**@file        DataSet.h
 * @brief       this file provide DataSet class
 * @details     DataSet class is a minimal dynamic container. This file includes
 *              - Node class as an abstract interface,
 *              - Shape class control the shape of the tensor
 *              - Tensor class inherited from Node
 *              - DataSet class with a tree structure for storage of Tensor
 *
 * @author      [author] [author2] [author3]
 * @date        [latest-date]
 * @version     [version]
 * @copyright   [copyright]
 **********************************************************************************
 * @par revision [logs]:
 * <table>
 * <tr><th> Date    <th> Version    <th> Author    <th> Description
 * <tr><td>[date]   <td>[version]   <td>[author]   <td> [commit]
 * </table>
 *
 **********************************************************************************
 */

#ifndef DataSet_H
#define DataSet_H

#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Exception.h"
#include "concat.h"
#include "types.h"

namespace PROJECT_NS {

/**
 * control the output printing format
 */
constexpr inline int FMT_WIDTH(int X) { return X + 7; }
#define FMT(X)                                                            \
    " " << std::setiosflags(std::ios::scientific) /*scientific notation*/ \
        << std::setprecision(X)                   /*precision*/           \
        << std::right                             /*alignment*/           \
        << std::setw(FMT_WIDTH(X))                /*width of text*/


/**********************************************************************************
    This file includes classes:
    1) Node
    2) Shape
    3) Tensor<T>
    4) DataSet
**********************************************************************************/

/**
 * Node class is an abstract interface. A Node will store a tensor or another Node
 */
class Node {
   public:
    using SizeType = std::size_t;
    using DataType = void;

    virtual ~Node() { _delete(); }  //< deconstructor

    virtual std::string repr() = 0;

    inline kids_dtype type() { return _type; }

    inline SizeType size() { return _size; }

    inline DataType* data() { return _data; }

   protected:
    friend class DataSet;

    kids_dtype _type = kids_void_type;
    SizeType _size   = 0;
    DataType* _data  = nullptr;
    bool _ownership  = false;

   private:
    virtual void _delete(){};  // virtual function called by deconstructor
};

/**
 * Shape class provide information about of a Tensor's shape
 */
class Shape {
   public:
    ///< construct from a vector
    Shape(std::vector<std::size_t> dims) : _rank{dims.size()}, _dims{dims}, _ldims(dims.size(), 0), _size{1} {
        _ldims[_ldims.size() - 1] = 1;
        _size                     = dims[_dims.size() - 1];
        for (int i = _dims.size() - 2; i >= 0; --i) {
            _ldims[i] = _ldims[i + 1] * (_dims[i + 1]);
            _size *= _dims[i];
        }
    }

    ///< construct from a number (rank-1 Shape)
    Shape(std::size_t size) : _rank{1}, _dims{{size}}, _ldims{1}, _size{size} {};

    /**
     * Delete new operators to prevent dynamic memory allocation
     */
    void* operator new(size_t)   = delete;
    void* operator new[](size_t) = delete;

    ///< get rank of a Shape
    inline int rank() const { return _rank; }

    ///< get data size described by a Shape
    inline int size() const { return _size; }

   private:
    std::size_t _rank;                ///< rank of the shape
    std::vector<std::size_t> _dims;   ///< dimensions for each rank
    std::vector<std::size_t> _ldims;  ///< leading dimensions
    std::size_t _size;                ///< size of data
};

/**
 * Tensor class is the container for array-like data
 */
template <typename T>
class Tensor final : public Node {
   public:
    using SizeType = std::size_t;
    using DataType = T;

    Tensor(Shape S) : _shape{S} {
        _type = as_enum<T>();
        _size = _shape.size();
        _data = (void*) new T[_size];
        memset(_data, 0, _size * sizeof(T));
        _ownership = true;
    }
    virtual ~Tensor() { _delete(); }

    virtual std::string repr() {
        std::ostringstream os;
        T* ptr = (T*) data();
        os << FMT(0) << as_str<T>();
        os << FMT(0) << _size;
        os << "\n";
        for (int i = 0; i < _size; ++i) os << FMT(8) << ptr[i];
        return os.str();
    }

   private:
    Shape _shape;

    virtual void _delete() {
        if (_ownership) delete[] static_cast<T*>(_data);
        _ownership = false;
    }
};

/**
 * DataSet class is a tree-structured container for storage of Tensor and
 * other DataSet.
 */
class DataSet final : public Node {
   private:
    DataSet& operator=(const DataSet&) = delete;

   public:
    using DataType = std::map<std::string, Node*>;

    DataSet() {
        _size      = 1;  // only one TensorBuffer will be newed
        _type      = kids_dataset_type;
        _data      = (void*) new DataType();
        _ownership = true;
    };

    virtual ~DataSet() {
        if (_ownership) _delete();
        _ownership = false;
    }

    template <typename T>
    T* def(const std::string& key, Shape S = 1) {
        DataSetKeyParser kh = DataSetKeyParser(key);
        DataType* d_ptr     = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) node = new DataSet;

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto& leaf_node = (*d_ptr)[kh.terms.back()];
        if (!leaf_node) leaf_node = new Tensor<T>(S);

        if (leaf_node->type() != as_enum<T>() || S.size() != leaf_node->size())
            throw std::runtime_error("doubly conflicted definition!");
        return static_cast<T*>(leaf_node->data());
    }

    template <typename T>
    T* def(const std::string& key, T* arr_in, Shape S = 1) {
        T* arr = def<T>(key, S);
        for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
        return arr;
    }

    template <typename T>
    T* def(const std::string& key, const std::string& key_in) {
        auto inode = node(key_in);
        if (inode->type() == kids_dataset_type)
            throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
        return def<T>(key, (T*) inode->data(), inode->size());
    }

    template <typename T>
    DataSet& _def(const std::string& key, Shape S = 1) {
        def<T>(key, S);
        return *this;
    }

    DataSet& _def(const std::string& key, const std::string& key_in) {
        auto leaf_node = node(key_in);
        switch (leaf_node->type()) {
            case kids_bool_type:
                def<kids_bool>(key, key_in);
                break;
            case kids_int_type:
                def<kids_int>(key, key_in);
                break;
            case kids_real_type:
                def<kids_real>(key, key_in);
                break;
            case kids_complex_type:
                def<kids_complex>(key, key_in);
                break;
            default:
                break;
        }
        return *this;
    }

    DataSet& _undef(const std::string& key) {
        DataSetKeyParser kh = DataSetKeyParser(key);
        DataType* d_ptr     = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) return *this;

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto it = d_ptr->find(kh.terms.back());
        if (it != d_ptr->end()) {
            delete it->second;
            d_ptr->erase(it);
        }
        return *this;
    }

    Node* node(const std::string& key) {
        DataSetKeyParser kh = DataSetKeyParser(key);
        DataType* d_ptr     = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) throw std::runtime_error(std::string{key} + " : access undefined key!");

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto& leaf_node = (*d_ptr)[kh.terms.back()];
        if (!leaf_node) throw std::runtime_error(std::string{key} + " : access undefined key!");
        return leaf_node;
    }

    template <typename T = DataSet>
    T* at(const std::string& key) {
        auto leaf_node = node(key);
        if (leaf_node->type() != as_enum<T>()) throw std::runtime_error("bad conversion!");
        if (leaf_node->type() == kids_dataset_type) return (T*) (leaf_node);
        return static_cast<T*>(leaf_node->data());
    }

    virtual std::string repr() {
        std::ostringstream os;
        std::vector<std::tuple<std::string, Node*>> stack;

        stack.push_back(std::make_tuple("", this));

        while (!stack.empty()) {
            auto [parent, currentNode] = stack.back();
            stack.pop_back();

            DataType* d_ptr = static_cast<DataType*>(currentNode->data());

            for (auto& i : (*d_ptr)) {
                std::string key = (parent == "") ? i.first : parent + "." + i.first;
                Node* inode     = i.second;
                if (inode->type() == kids_dataset_type) {
                    stack.push_back(std::make_tuple(key, inode));
                } else {
                    os << key << "\n" << inode->repr() << "\n\n";
                }
            }
        }
        return os.str();
    }

    virtual void dump(std::ostream& os) { os << repr(); }

    virtual void load(std::istream& is) {
        std::string key, typeflag;
        int size;
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

   private:
    class DataSetKeyParser {
       public:
        std::vector<std::string> terms;
        DataSetKeyParser(const std::string& key, const std::string& delimiter = ".") {
            size_t start = 0, end;
            while ((end = key.find(delimiter, start)) != std::string::npos) {
                terms.emplace_back(key, start, end - start);
                start = end + delimiter.length();
            }
            terms.emplace_back(key, start);
        }
    };

    virtual void _delete() {
        DataType* d_ptr = static_cast<DataType*>(_data);
        for (auto& i : (*d_ptr)) delete i.second;
        delete d_ptr;  // essential!
    }
};

};  // namespace PROJECT_NS

#endif  // DataSet_H

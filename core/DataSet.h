/**@file        DataSet.h
 * @brief       Declaration of the DataSet class and related classes.
 * @details     This file provides the declaration of the DataSet class. It serves
 *              as a minimal dynamic container for storing tensors. It also
 *              includes declarations for the following supporting classes:
 *              - Node: An abstract interface class for tree nodes.
 *              - Tensor: A class inherited from Node, representing tensors.
 *              - DataSet: A class implementing a tree structure for the storage
 *
 * @author      Xin He
 * @date        2024-03
 * @version     1.0
 * @copyright   GNU Lesser General Public License (LGPL)
 *
 *              Copyright (c) 2024 Xin He, Liu-Group
 *
 *  This software is part of the research conducted by the Prof. Liu's Group at the
 *  College of Chemistry and Molecular Engineering (CCME), Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
 * @todo
 *  Implement a smart pointer version for the DataSet class. One possible approach
 *  is to use `std::unique_ptr<T[]>`, as demonstrated in the file:
 *      `DataSet.cppversion3_smartlymanage`.
 *  However, it should be noted that accessing data is not currently handled
 *  correctly with this implementation. On the other hand, `std::unique_ptr` is
 *  unsuitable for scenarios requiring shared access to the storage of the DataSet
 *  object. While `std::shared_ptr<T[]>` could potentially address this issue, it
 *  is not as well-supported in C++. Consequently, designing the DataSet class with
 *  `std::shared_ptr` is necessary.
 *
 *  Due to this limitation, storage is managed manually in current version. Copying
 *  a DataSet to another instance is prohibited, and reassigning a DataSet is also
 *  disallowed.  Additionally, there is currently no interface provided for Python
 *  for reassign a DataSet object.
 *
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-03-29  <td> initial version. Seperate Shape class to another file.
 *                          Add help() function. Review more over smart pointers.
 *                          Improve print format.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_DataSet_H
#define KIDS_DataSet_H

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
#include "Shape.h"
#include "Types.h"
#include "concat.h"

namespace PROJECT_NS {

/**
 * control the io printing format
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
    2) Tensor<T>
    3) DataSet
**********************************************************************************/

/**
 * Node class is an abstract interface. A Node will store a tensor or other Nodes
 */
class Node {
   public:
    using SizeType = std::size_t;
    using DataType = void;

    virtual ~Node() { _delete(); }  //< deconstructor

    virtual std::string repr() = 0;

    virtual std::string help(const std::string& name) = 0;

    inline kids_dtype type() { return _type; }

    inline SizeType size() { return _size; }

    inline DataType* data() { return _data; }

   protected:
    friend class DataSet;

    kids_dtype _type      = kids_void_type;
    SizeType   _size      = 0;
    DataType*  _data      = nullptr;
    bool       _ownership = false;

   private:
    virtual void _delete(){};  // virtual function called by deconstructor
};

/**
 * Tensor class is the container for array-like data
 */
template <typename T>
class Tensor final : public Node {
   public:
    using SizeType = std::size_t;
    using DataType = T;

    Tensor(Shape S, const std::string& info = "") : _shape{S}, _doc_info{info} {
        _type = as_enum<T>();
        _size = _shape.size();
        _data = (void*) new T[_size];
        memset(_data, 0, _size * sizeof(T));
        _ownership = true;
    }
    virtual ~Tensor() { _delete(); }

    virtual std::string repr() {
        std::ostringstream os;
        T*                 ptr = (T*) data();
        os << as_str<T>();
        os << FMT(0) << _size;
        os << "\n";
        for (int i = 0; i < _size; ++i) os << FMT(8) << ptr[i];
        return os.str();
    }

    virtual std::string help(const std::string& name) { return _doc_info; }

   private:
    Shape       _shape;
    std::string _doc_info;

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
    T* def(const std::string& key, Shape S = 1, const std::string& info = "") {
        DataSetKeyParser kh    = DataSetKeyParser(key);
        DataType*        d_ptr = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) node = new DataSet;

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto& leaf_node = (*d_ptr)[kh.terms.back()];
        if (!leaf_node) leaf_node = new Tensor<T>(S, info);

        if (leaf_node->type() != as_enum<T>() || S.size() != leaf_node->size())
            throw std::runtime_error("doubly conflicted definition!");
        return static_cast<T*>(leaf_node->data());
    }

    template <typename T>
    T* def(const std::string& key, T* arr_in, Shape S = 1, const std::string& info = "") {
        T* arr = def<T>(key, S, info);
        for (int i = 0; i < S.size(); ++i) arr[i] = arr_in[i];
        return arr;
    }

    template <typename T>
    T* def(const std::string& key, const std::string& key_in, const std::string& info = "") {
        auto inode = node(key_in);
        if (inode->type() == kids_dataset_type)
            throw std::runtime_error(std::string{key_in} + " : failed copying dataset");
        return def<T>(key, (T*) inode->data(), inode->size(), info);
    }

    template <typename T>
    DataSet& _def(const std::string& key, Shape S = 1, const std::string& info = "") {
        def<T>(key, S, info);
        return *this;
    }

    DataSet& _def(const std::string& key, const std::string& key_in, const std::string& info = "") {
        auto leaf_node = node(key_in);
        switch (leaf_node->type()) {
            case kids_bool_type:
                def<kids_bool>(key, key_in, info);
                break;
            case kids_int_type:
                def<kids_int>(key, key_in, info);
                break;
            case kids_real_type:
                def<kids_real>(key, key_in, info);
                break;
            case kids_complex_type:
                def<kids_complex>(key, key_in, info);
                break;
            default:
                break;
        }
        return *this;
    }

    DataSet& _undef(const std::string& key) {
        DataSetKeyParser kh    = DataSetKeyParser(key);
        DataType*        d_ptr = static_cast<DataType*>(_data);

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
        DataSetKeyParser kh    = DataSetKeyParser(key);
        DataType*        d_ptr = static_cast<DataType*>(_data);

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

    virtual std::string help(const std::string& name) {
        std::ostringstream                          os;
        std::vector<std::tuple<std::string, Node*>> stack;

        Node* inode = this;
        if (name != "") inode = node(name);

        stack.push_back(std::make_tuple("", inode));
        while (!stack.empty()) {
            auto [parent, currentNode] = stack.back();
            stack.pop_back();

            DataType* d_ptr = static_cast<DataType*>(currentNode->data());

            for (auto& i : (*d_ptr)) {
                std::string key   = (parent == "") ? i.first : parent + "." + i.first;
                Node*       inode = i.second;
                if (inode->type() == kids_dataset_type) {
                    stack.push_back(std::make_tuple(key, inode));
                } else {
                    os << key << ":\n\t" << inode->help("") << "\n";
                }
            }
        }
        return os.str();
    }

    virtual std::string repr() {
        std::ostringstream                          os;
        std::vector<std::tuple<std::string, Node*>> stack;

        stack.push_back(std::make_tuple("", this));

        while (!stack.empty()) {
            auto [parent, currentNode] = stack.back();
            stack.pop_back();

            DataType* d_ptr = static_cast<DataType*>(currentNode->data());

            for (auto& i : (*d_ptr)) {
                std::string key   = (parent == "") ? i.first : parent + "." + i.first;
                Node*       inode = i.second;
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

#endif  // KIDS_DataSet_H

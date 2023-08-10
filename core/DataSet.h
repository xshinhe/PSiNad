/**
 * @file DataSet.h
 * @author xshinhe
 * @date 2023-05
 * @brief class for `DataSet` object (alias of Tree)
 * @details
 *  It stores variables in a variable-list, with almost zero cost in copy
 *  (i.e. fetch the pointer of the data only). It support a probability that
 *  different kernels/solvers could access the same data storage.
 *  It also facilitates to save current status into a file, and fully reload
 *  it agian whenever you want.
 * @structure
 *  this file realizes three class.
 *  1) NodeGeneric: interface can contain any type of data by (void*)
 *  2) Node<T>: leaf NodeGeneric, which contains an array of type T.
 *  3) Tree: non-leaf NodeGeneric, which contains `Tree::NodeBuffer` type.
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

// format control for dumping the state
inline int FMT_WIDTH_SIZE(int X) { return X + 9; }
#define FMT(X) \
    " " << std::setiosflags(std::ios::scientific) << std::setprecision(X) << std::right << std::setw(FMT_WIDTH_SIZE(X))

namespace PROJECT_NS {

template <typename T>
struct DSVar {
    T* ptr;
    std::string key;
    int size;
};

namespace details {  // hidden scope

class Tree;  // declaration

/**
 * @brief Generic Node interface
 */
class NodeGeneric {
   public:
    enum class Format {
        Free,  // free format
        Json   // json format
    };
    enum class NodeType {
        Void,     // nothing
        Int,      // int c-array
        Real,     // num_real c-array
        Complex,  // num_complex c-array
        Tree      // tree
    };
    using SizeType = std::size_t;
    using DataType = void;
    using InfoType = std::tuple<NodeType, SizeType, DataType*, NodeGeneric*>;
    using HelpType = std::tuple<NodeType, std::string>;

    template <typename T>
    static inline HelpType TypeHelper();

    virtual ~NodeGeneric() { _delete(); }  //< deconstructor

    virtual std::string repr(Format format, const std::string& lead = "") = 0;

    inline NodeType type() { return _type; }

    inline SizeType size() { return _size; }

    inline DataType* data() { return _data; }

    inline InfoType info() { return std::make_tuple(_type, _size, _data, this); }

   protected:
    NodeType _type  = NodeType::Void;
    SizeType _size  = 0;
    DataType* _data = nullptr;

   private:
    virtual void _delete(){};  // virtual function called by deconstructor
};

#define DEFINE_TYPE_BIND_WITH_NODETYPE(T, NT)                   \
    template <>                                                 \
    inline NodeGeneric::HelpType NodeGeneric::TypeHelper<T>() { \
        return std::make_tuple(NodeGeneric::NodeType::NT, #NT); \
    }

DEFINE_TYPE_BIND_WITH_NODETYPE(int, Int);
DEFINE_TYPE_BIND_WITH_NODETYPE(num_real, Real);
DEFINE_TYPE_BIND_WITH_NODETYPE(num_complex, Complex);
DEFINE_TYPE_BIND_WITH_NODETYPE(Tree, Tree);

/**
 * @brief leaf NodeGeneric, which contains an array of type T
 * @tparam T customized type for Node (in int, num_real, num_complex)
 */
template <typename T>
class Node final : public NodeGeneric {
   public:
    using SizeType = std::size_t;
    using DataType = T;

    Node(const std::size_t& size) {
        _type = std::get<0>(NodeGeneric::TypeHelper<T>());
        _size = size;
        _data = (void*) new T[_size];
        memset(_data, 0, _size * sizeof(T));
        _init = true;
    }
    virtual ~Node() { _delete(); }

    virtual std::string repr(NodeGeneric::Format format, const std::string& lead = "") {
        std::ostringstream os;
        T* ptr = (T*) data();
        switch (format) {
            case NodeGeneric::Format::Free: {
                os << FMT(0) << std::get<1>(NodeGeneric::TypeHelper<T>());
                os << FMT(0) << _size;
                os << FMT(8) << ptr[0];
                for (int i = 1; i < _size; ++i) os << FMT(8) << ptr[i];
                break;
            }
            case NodeGeneric::Format::Json: {
                os << "[";
                os << FMT(8) << ptr[0];
                for (int i = 1; i < _size; ++i) os << "," << FMT(8) << ptr[i];
                os << "]";
                break;
            }
        }
        return os.str();
    }

   private:
    friend class Tree;
    bool _init = false;

    virtual void _delete() {
        if (_init) delete[] static_cast<T*>(_data);
        _init = false;
    }
};

/**
 * @brief non-leaf NodeGeneric, which contains `Tree::NodeBuffer` type.
 */
class Tree : public NodeGeneric {
   private:
    Tree& operator=(const Tree&) = delete;

   public:
    using KeyNode     = std::tuple<std::string, NodeGeneric*>;
    using KeyNodeList = std::vector<KeyNode>;
    using KeyNodeMap  = std::map<std::string, NodeGeneric*>;
    using DataType    = KeyNodeMap;

    using Generic = NodeGeneric;
    using Format  = NodeGeneric::Format;
    using Type    = NodeGeneric::NodeType;

    /**
     * @brief constructor: to allocate an instance of NodeBuffer (size=1)
     */
    Tree() {
        _size     = 1;  // only one NodeBuffer will be newed
        _type     = NodeType::Tree;
        _data     = (void*) new DataType();
        ownership = true;
    };

    Tree(const Tree& itree) {     // always referenced-copy (deep copy see at: shapelike())
        _size     = itree._size;  // only one NodeBuffer will be newed
        _type     = itree._type;
        _data     = itree._data;
        ownership = false;
    }

    /**
     * @brief deconstructor
     */
    virtual ~Tree() {
        if (ownership) _delete();
    }

    /**
     * @brief      define a leaf-node of in tree storing c-array
     *
     * @param[in]  key   The key, i.e., the address of the c-array, adjointed by dot delimiter
     * @param[in]  size  The size of the array
     *
     * @tparam     T     The type of the data (int, num_real, num_complex)
     *
     * @return     Tree reference after add definition
     */
    template <typename T>
    Tree& def(const std::string& key, std::size_t size = 1) {
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);
        if (kh.is_leaf) {  // leaf node
            if (d_ptr->find(key) != d_ptr->end()) throw state_conflicted_key_error(key);
            (*d_ptr)[key] = new Node<T>(size);
        } else {  // non-leaf node
            if (d_ptr->find(kh.key1) == d_ptr->end()) (*d_ptr)[kh.key1] = new Tree;
            if ((*d_ptr)[kh.key1]->type() != NodeType::Tree) throw state_conflicted_key_error(kh.key1);
            ((Tree*) (*d_ptr)[kh.key1])->def<T>(kh.key2, size);
        }
        return *this;
    }

    /**
     * @brief      undefine leaf/non-leaf node(s) in a tree
     *
     * @param[in]  key   The key
     *
     * @return     Tree reference after undefintion
     */
    Tree& undef(const std::string& key) {
        KeyHelper kh    = KeyHelper(key, false);
        Tree* psubtree  = kh.is_leaf ? this : &subref(kh.key1);
        DataType* d_ptr = static_cast<DataType*>(psubtree->_data);
        if (d_ptr->find(kh.key2) == d_ptr->end()) throw state_undefined_key_error(kh.key2);
        delete (*d_ptr)[kh.key2];
        d_ptr->erase(kh.key2);
        return *this;
    }


    /**
     * @brief      get a sub-tree of a tree
     *
     * @param[in]  key   The key of a subtree
     *
     * @return     Tree reference
     */
    Tree& subref(const std::string& key) {  // @bugs it will destroy the alias structure
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);
        if (kh.is_leaf) {
            if (d_ptr->find(key) == d_ptr->end()) throw state_undefined_key_error(key);
            if (NodeType::Tree != (*d_ptr)[key]->type()) throw state_mismatched_type_error(key);
            return *((Tree*) (*d_ptr)[key]);
        }
        if (d_ptr->find(kh.key1) == d_ptr->end()) throw state_undefined_key_error(kh.key1);
        return ((Tree*) (*d_ptr)[kh.key1])->subref(kh.key2);
    };

    /**
     * @brief get information (type & size & storage) of a key in the tree
     */
    NodeGeneric::InfoType info(const std::string& key) noexcept {
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);
        if (kh.is_leaf && d_ptr->find(key) != d_ptr->end()) {
            NodeGeneric* inode = (*d_ptr)[key];
            return inode->info();
        } else if (!kh.is_leaf && d_ptr->find(kh.key1) != d_ptr->end()) {
            return ((Tree*) (*d_ptr)[kh.key1])->info(kh.key2);
        }
        return std::make_tuple(NodeGeneric::NodeType::Void, 0, nullptr, nullptr);
    };

    /**
     * @brief ask if the key exists in the tree structure
     */
    bool has_key(const std::string& key) { return (std::get<3>(info(key)) != nullptr); }

    /**
     * @brief register an array and return its pointer

     * @param key : the address of the array, adjointed by dot delimiter
     * @param size_array : the size of the array
     * @tparam T : the type of the data
     * @return T* pointer.
     *     1) If the key of array has been existed, it will only return the
     *     pointer of the array;
     *     2) otherwise, it will define a new array assoicated with the key,
     *     and further return its pointer.
     */
    template <typename T>
    T* reg(const std::string& key, std::size_t size_array = 1, bool required = false) {
        if (!has_key(key)) def<T>(key, size_array);
        auto&& res = info(key);
        if (std::get<0>(NodeGeneric::TypeHelper<T>()) != std::get<0>(res)) throw state_mismatched_type_error(key);
        if (size_array > 0 && size_array != std::get<1>(res)) throw state_mismatched_size_error(key);
        return (T*) std::get<2>(res);
    };

    template <typename T>
    void set(const std::string& key, T* array, std::size_t size_array = 1) {
        if (!has_key(key)) def<T>(key, size_array);
        auto&& res = info(key);
        if (std::get<0>(NodeGeneric::TypeHelper<T>()) != std::get<0>(res)) throw state_mismatched_type_error(key);
        if (size_array > 0 && size_array != std::get<1>(res)) throw state_mismatched_size_error(key);
        T* ptr = (T*) std::get<2>(res);
        for (int i = 0; i < size_array; ++i) ptr[i] = array[i];
    };

    KeyNodeList flatten(bool indeep = false, const std::string& parent = "") {
        KeyNodeList list;
        DataType* d_ptr = static_cast<DataType*>(_data);
        for (auto& i : (*d_ptr)) {
            std::string key    = (parent == "") ? i.first : utils::concat(parent, ".", i.first);
            NodeGeneric* inode = i.second;
            if (inode->type() == NodeGeneric::NodeType::Tree) {
                auto&& list2 = ((Tree*) inode)->flatten(indeep, key);
                list.insert(list.end(), std::make_move_iterator(list2.begin()), std::make_move_iterator(list2.end()));
            } else {
                list.push_back(std::make_tuple(key, inode));
            }
        }
        return list;
    }

    virtual std::string repr(NodeGeneric::Format format = NodeGeneric::Format::Free, const std::string& lead = "") {
        std::ostringstream os;
        DataType* d_ptr = static_cast<DataType*>(_data);
        switch (format) {
            case NodeGeneric::Format::Free: {
                for (auto& i : flatten()) { os << std::get<0>(i) << std::get<1>(i)->repr(format, lead) << "\n"; }
                break;
            }
            case NodeGeneric::Format::Json: {
                os << lead << "{\n";
                std::string nextlead = lead + "  ";
                for (auto& i : (*d_ptr)) {
                    NodeGeneric* inode = i.second;
                    os << nextlead << "\"" << i.first << "\" : ";
                    os << nextlead << inode->repr(format, nextlead);
                    os << ",\n";
                }
                os << lead << "}";
                break;
            }
        }
        return os.str();
    }

    /**
     * @brief dump Tree information to a filestream or a file
     */
    virtual void dump(std::ostream& os, NodeGeneric::Format format = NodeGeneric::Format::Free) { os << repr(format); }

    virtual void dump(const std::string& file, NodeGeneric::Format format = NodeGeneric::Format::Free) {
        std::ofstream ofs(file);
        dump(ofs, format);
        ofs.close();
    }

    /**
     * @brief load Tree information from a filestream or a file
     */
    virtual void load(std::istream& is, NodeGeneric::Format format = NodeGeneric::Format::Free) {
        switch (format) {
            case NodeGeneric::Format::Free: {  // as an easily implemented format
                std::string key, typeflag;
                int size;
                while (is >> key >> typeflag >> size) {
                    if (typeflag == "Int") {
                        int* ptr = reg<int>(key, size);
                        for (int i = 0; i < size; ++i) is >> ptr[i];
                    } else if (typeflag == "Real") {
                        num_real* ptr = reg<num_real>(key, size);
                        for (int i = 0; i < size; ++i) is >> ptr[i];
                    } else if (typeflag == "Complex") {
                        num_complex* ptr = reg<num_complex>(key, size);
                        for (int i = 0; i < size; ++i) is >> ptr[i];
                    }
                }
                break;
            }
            case NodeGeneric::Format::Json: {
                throw basic_error("Not implemented");
                break;
            }
        }
    }

    virtual void load(const std::string& file, NodeGeneric::Format format = NodeGeneric::Format::Free) {
        std::ifstream ifs(file);
        load(ifs, format);
        ifs.close();
    }

   private:
    bool ownership = true;

    /**
     * @brief      This Class parses key sequences.
     */
    class KeyHelper {
       public:
        bool is_leaf;
        std::string key1, key2;
        KeyHelper(const std::string key, bool from_beginning = true) {
            std::string::size_type ipos = from_beginning ? key.find_first_of('.') : key.find_last_of('.');
            is_leaf                     = (ipos == std::string::npos);
            if (is_leaf) {
                key1 = from_beginning ? key : "";
                key2 = from_beginning ? "" : key;
            } else {  // find a delimiter
                key1 = key.substr(0, ipos++);
                key2 = key.substr(ipos, key.length() - ipos);
            }
        }
    };

    /**
     * @brief implemention of deconstructor
     */
    virtual void _delete() {
        DataType* d_ptr = static_cast<DataType*>(_data);
        for (auto& i : (*d_ptr)) delete i.second;
        delete d_ptr;  // essential!
    }
};

};  // namespace details

using DataSet = details::Tree;  // exported components

};  // namespace PROJECT_NS

#endif  // DataSet_H

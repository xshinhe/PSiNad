#ifndef DataSetTensor_H
#define DataSetTensor_H

#include <algorithm>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

namespace PROJECT_NS {

using kids_real    = double;
using kids_complex = std::complex<double>;

// format control for dumping the state
constexpr inline int FMT_WIDTH(int X) { return X + 7; }
#define FMT(X) \
    " " << std::setiosflags(std::ios::scientific) << std::setprecision(X) << std::right << std::setw(FMT_WIDTH(X))

/////////////////////////////////////////////////////////
// forward declaration of intrinsic class
class Dimen;          ///< dynamics Dimenson class
class Shape;          ///< dynamics Shape class
class VARIABLE_BASE;  ///< interface for variables
template <class T>    ///
class VARIABLE;       ///< wrapper for variables
class Node;           ///< interface for storage
template <class T>    ///
class Tensor;         ///< storage of tensor
class DataSet;        ///< storage of a dataset
/////////////////////////////////////////////////////////

/**
 * @brief      This class describes a dynamic dimension.
 * @usage
 *      1) `Dimen dim;` define a dim
 *      2) `dim.set(10);` assignment the value (not completed)
 *      3) `dim = 10;` assignment the value (completed and can be only called once)
 *      4) `dim.is_completed();` check if dimension is completed
 *      5) `dim();` get the value
 */
class Dimen {
   public:
    Dimen() : _val{0}, _completed{false} {};
    Dimen(std::size_t val) : _val{val}, _completed{true} {};

    // Delete new operators to prevent dynamic memory allocation
    void* operator new(std::size_t size)   = delete;
    void* operator new[](std::size_t size) = delete;

    inline int operator()() const { return _val; }

    inline bool is_completed() const { return _completed; }

    inline void set(int val, bool completed = false);

    void operator=(int val) { set(val, true); }

   private:
    friend class Shape;

    bool _completed;
    std::size_t _val;
    std::vector<Shape*> associated_shapes;
};

/**
 * @brief      This class describes a tensor shape from the dynamic dimension
 * @usage
 *      1) `Shape S12({&dim1, &dim2});`
 *      2) `S12();` total size
 *      3) `S12(0)`; size along axis-0
 *      4) `dim1.set(newvalue);` dynamic update of Shape S12
 */
class Shape {
   public:
    Shape(std::vector<Dimen*> dims)
        : nrank{dims.size()},    // Rank of the shape
          _dims{dims},           // Dimensions for each rank
          _ldas(dims.size(), 0)  // Leading dimensions for each rank
    {
        // Associate this shape with each dimension
        for (auto& dim : _dims) dim->associated_shapes.push_back(this);
        _update();
    }

    // Delete new operators to prevent dynamic memory allocation
    void* operator new(size_t size)   = delete;
    void* operator new[](size_t size) = delete;

    inline int operator()() const { return _totalsize; }

    inline int operator()(std::size_t idxr) const { return _dims[idxr]->_val; }

    inline bool is_completed() const { return _completed; }

    inline int rank() const { return nrank; }

   private:
    friend class Dimen;

    bool _completed;

    std::vector<Dimen*> _dims;
    std::vector<std::size_t> _ldas;
    std::size_t nrank, _totalsize;

    // Update dimensions and leading dimensions
    void _update() {
        if (std::any_of(_dims.begin(), _dims.end(), [](const Dimen* idim) { return !idim->_completed; })) {
            _completed = false;
            return;  // Do nothing unless all Dimen objects are completed
        }

        // Calculate leading dimensions and total size
        _ldas[_ldas.size() - 1] = 1;
        for (int i = _dims.size() - 2; i >= 0; --i) _ldas[i] = _ldas[i + 1] * _dims[i]->_val;
        _totalsize = _dims[0]->_val * _ldas[0];
        _completed = true;
    }
};

// appending realization: Assign value to Dimen object
void Dimen::set(int val, bool completed) {
    if (_completed) return;
    _val       = val;
    _completed = completed;
    for (auto& ishape : associated_shapes) ishape->_update();
}

class VARIABLE_BASE {
   public:
    virtual std::string name() = 0;
    virtual std::string help() = 0;

    static std::vector<VARIABLE_BASE*> _LIST;
};
std::vector<VARIABLE_BASE*> VARIABLE_BASE::_LIST;

template <class T>
class VARIABLE final : public VARIABLE_BASE {
   public:
    VARIABLE(const std::string& name, Shape* shape, const std::string& doc) : _name{name}, _shape{shape}, _doc{doc} {
        VARIABLE_BASE::_LIST.push_back(this);
    }

    std::string help() { return _doc; }
    std::string name() { return _name; }
    std::size_t size() { return (*_shape)(); }
    Shape* shape() const { return _shape; }

   private:
    std::string _name;
    T* _data;
    Shape* _shape;
    std::string _doc;
    bool allocated = false;
};

/**
 * @brief      This class describes a node (the interface).
 */
class Node {
   public:
    enum class NodeType { Void, Bool, Int, Real, Complex, DataSet };
    using SizeType = std::size_t;
    using DataType = void;

    virtual ~Node() { _delete(); }  //< deconstructor

    virtual std::string repr() = 0;

    inline NodeType type() { return _type; }

    inline SizeType size() { return _size; }

    inline DataType* data() { return _data; }

   protected:
    NodeType _type  = NodeType::Void;
    SizeType _size  = 0;
    DataType* _data = nullptr;
    bool _ownership = false;

   private:
    virtual void _delete(){};  // virtual function called by deconstructor
};

///< binding for types with enums and strings
template <typename T>
Node::NodeType AsNodeType();
template <typename T>
std::string AsString();
#define DEFINE_TYPE_BIND_WITH_NODETYPE(T, NT) \
    template <>                               \
    inline Node::NodeType AsNodeType<T>() {   \
        return Node::NodeType::NT;            \
    }                                         \
    template <>                               \
    inline std::string AsString<T>() {        \
        return #NT;                           \
    }
DEFINE_TYPE_BIND_WITH_NODETYPE(int, Int);
DEFINE_TYPE_BIND_WITH_NODETYPE(bool, Bool);
DEFINE_TYPE_BIND_WITH_NODETYPE(kids_real, Real);
DEFINE_TYPE_BIND_WITH_NODETYPE(kids_complex, Complex);
DEFINE_TYPE_BIND_WITH_NODETYPE(DataSet, DataSet);


/**
 * @brief      This class describes a tensor.
 *
 * @tparam     T     { value type }
 */
template <typename T>
class Tensor final : public Node {
   public:
    using SizeType = std::size_t;
    using DataType = T;

    Tensor(const std::size_t& size) {
        _type  = AsNodeType<T>();
        _size  = size;
        _data  = (void*) new T[_size];
        _shape = nullptr;
        memset(_data, 0, _size * sizeof(T));
        _ownership = true;
    }
    Tensor(Shape* shape_ptr) {
        if (!shape_ptr->is_completed()) throw std::runtime_error("uncompleted shape");

        _type  = AsNodeType<T>();
        _size  = (*shape_ptr)();
        _data  = (void*) new T[_size];
        _shape = shape_ptr;
        memset(_data, 0, _size * sizeof(T));
        _ownership = true;
    }
    // Tensor(VARIABLE<T>& VAR) : Tensor(VAR.shape()){};

    virtual ~Tensor() { _delete(); }

    virtual std::string repr() {
        std::ostringstream os;
        T* ptr = (T*) data();
        os << FMT(0) << AsString<T>();
        os << FMT(0) << _size;
        if (_shape == nullptr) {
            os << FMT(0) << 1 << FMT(0) << _size;
        } else {
            os << FMT(0) << (_shape)->rank();
            for (int k = 0; k < _shape->rank(); ++k) os << FMT(0) << (*_shape)(k);
        }
        os << std::endl;
        for (int i = 0; i < _size; ++i) os << FMT(8) << ptr[i];
        return os.str();
    }

   private:
    friend class DataSet;

    Shape* _shape;

    virtual void _delete() {
        if (_ownership) delete[] static_cast<T*>(_data);
        _ownership = false;
    }
};

inline int wrapsize(int i) { return i; }
inline int wrapsize(Shape* S) { return (*S)(); }

class DataSet final : public Node {
   private:
    DataSet& operator=(const DataSet&) = delete;

   public:
    using DataType = std::map<std::string, Node*>;

    DataSet() {
        _size      = 1;  // only one TensorBuffer will be newed
        _type      = Node::NodeType::DataSet;
        _data      = (void*) new DataType();
        _ownership = true;
    };

    virtual ~DataSet() {
        if (_ownership) _delete();
        _ownership = false;
    }

    template <typename T, typename U>
    T* def(const std::string& key, U TensorShape) {
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) node = new DataSet;

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto& leaf_node = (*d_ptr)[kh.terms.back()];
        if (!leaf_node) leaf_node = new Tensor<T>(TensorShape);

        if (leaf_node->type() != AsNodeType<T>() || wrapsize(TensorShape) != leaf_node->size())
            throw std::runtime_error("doubly conflicted definition!");
        return static_cast<T*>(leaf_node->data());
    }

    template <typename T>
    T* def(VARIABLE<T>& VAR) {
        return def<T>(VAR.name(), VAR.shape());
    }

    template <typename T = DataSet>
    T* at(const std::string& key) {
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) throw std::runtime_error("access undefined key!");

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto& leaf_node = (*d_ptr)[kh.terms.back()];
        if (!leaf_node) throw std::runtime_error("access undefined key!");

        if (leaf_node->type() != AsNodeType<T>()) throw std::runtime_error("bad conversion!");
        if (leaf_node->type() == Node::NodeType::DataSet) return (T*) (leaf_node);
        return static_cast<T*>(leaf_node->data());
    }

    void undef(const std::string& key) {
        KeyHelper kh    = KeyHelper(key);
        DataType* d_ptr = static_cast<DataType*>(_data);

        DataSet* currentNode = this;
        for (size_t i = 0; i < kh.terms.size() - 1; ++i) {
            auto& node = (*d_ptr)[kh.terms[i]];
            if (!node) return;

            currentNode = static_cast<DataSet*>(node);
            d_ptr       = static_cast<DataType*>(currentNode->_data);
        }

        auto it = d_ptr->find(kh.terms.back());
        if (it != d_ptr->end()) {
            delete it->second;
            d_ptr->erase(it);
        }
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
                std::string key = (parent == "") ? i.first : parent + "::" + i.first;
                Node* inode     = i.second;
                if (inode->type() == Node::NodeType::DataSet) {
                    stack.push_back(std::make_tuple(key, inode));
                } else {
                    os << key << inode->repr() << "\n";
                }
            }
        }
        return os.str();
    }

    virtual void dump(std::ostream& os) { os << repr(); }

    virtual void load(std::istream& is) {
        std::string key, typeflag;
        int size, rank, tmp;
        while (is >> key >> typeflag >> size) {
            is >> rank;
            for (int i = 0; i < rank; ++i) is >> tmp;
            if (typeflag == AsString<int>()) {
                int* ptr = def<int>(key, size);
                for (int i = 0; i < size; ++i) is >> ptr[i];
            } else if (typeflag == AsString<kids_real>()) {
                kids_real* ptr = def<kids_real>(key, size);
                for (int i = 0; i < size; ++i) is >> ptr[i];
            } else if (typeflag == AsString<kids_complex>()) {
                kids_complex* ptr = def<kids_complex>(key, size);
                for (int i = 0; i < size; ++i) is >> ptr[i];
            }
        }
    }

   private:
    class KeyHelper {
       public:
        std::vector<std::string> terms;
        KeyHelper(const std::string& key, const std::string& delimiter = "::") {
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

#define DATASET_DECLARE_VARIABLE(type, name) \
    namespace name {                         \
    VARIABLE<type> var;                      \
    };

#define DATASET_DEFINE_VARIABLE(type, name, shape, doc) \
    namespace name {                                    \
    VARIABLE<type> var(#name, shape, doc);              \
    };

};  // namespace PROJECT_NS

using namespace PROJECT_NS;

Dimen dim1;
Dimen dim2;
Shape S12({&dim1, &dim2});

DATASET_DEFINE_VARIABLE(double, integrator::x, &S12, "it is the coordinate");
DATASET_DEFINE_VARIABLE(double, integrator::p, &S12, "it is the momentum");
DATASET_DEFINE_VARIABLE(kids_complex, integrator::c, &S12, "it is the amplititude");

int main() {
    dim1 = 10;
    dim2 = 20;
    DataSet DS;
    DS.def<int>("a::b", 10);
    DS.def<double>("a::c::1", &S12);
    DS.def<double>("a::c::d", &S12);
    std::cout << DS.repr() << "\n";
    std::cout << DS.at("a::c")->repr() << "\n";
    std::cout << *DS.at<double>("a::c::1") << "\n";
    DS.undef("a::c");

    DS.def(integrator::x::var);
    DS.def(integrator::p::var);

    std::cout << DS.repr() << "\n";

    std::ofstream ofs{"test.ds"};
    DS.dump(ofs);
    ofs.close();

    std::ifstream ifs{"test.ds"};
    DataSet DS2;
    DS2.load(ifs);
    ifs.close();
    std::cout << DS2.repr() << "\n";

    for (auto& i : PROJECT_NS::VARIABLE_BASE::_LIST) { std::cout << i->name() << " : " << i->help() << "\n"; }
    return 0;
}

#endif  // DataSetTensor_H

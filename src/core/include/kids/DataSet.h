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
 *  This software is a product of Xin's PhD research conducted by Professor Liu's
 *  Group at the College of Chemistry and Molecular Engineering, Peking University.
 *  All rights are reserved by Peking University.
 *  You should have received a copy of the GNU Lesser General Public License along
 *  with this software. If not, see <https://www.gnu.org/licenses/lgpl-3.0.en.html>
 **********************************************************************************
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
 * <tr><td> 2024-04-18  <td> move to smart pointers. Instead of shared_ptr<T[]>, we
 *                          use shared_ptr<vector<T>> to manage memory.
 *                          However, there are limits:
 *                          1) vector<bool> is not STL, so it's not supported as a
 *                              part of Tensor<T>.
 *                          2) if vector is dynamically extended, the head pointer
 *                              may be changed.
 * <tr><td> 2024-04-23  <td> move to smart pointers. Using shared_ptr<T> and custom
 *                              delete function to help to manage an array of data.
 *
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_DataSet_H
#define KIDS_DataSet_H

#include <fstream>
#include <map>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

#include "kids/Exception.h"
#include "kids/Node.h"
#include "kids/Shape.h"
#include "kids/Tensor.h"
#include "kids/Types.h"
#include "kids/Variable.h"
#include "kids/concat.h"

namespace PROJECT_NS {

/**
 * DataSet class is a tree-structured container for storage of Tensor and
 * other DataSet.
 */
class DataSet final : public Node {
   private:
    DataSet& operator=(const DataSet&) = delete;

   public:
    using DataType = std::map<std::string, std::shared_ptr<Node>>;
    std::shared_ptr<DataType> _data;

    /**
     * Constructor for DataSet.
     */
    DataSet();

    /**
     * Define a variable of type kids_int.
     * @param var The variable to define.
     * @return Pointer to the defined variable.
     */
    kids_int* def(VARIABLE<kids_int>& var);

    /**
     * Define a variable of type kids_real.
     * @param var The variable to define.
     * @return Pointer to the defined variable.
     */
    kids_real* def(VARIABLE<kids_real>& var);

    /**
     * Define a variable of type kids_complex.
     * @param var The variable to define.
     * @return Pointer to the defined variable.
     */
    kids_complex* def(VARIABLE<kids_complex>& var);

    /**
     * Define an integer variable with a specified key, shape, and info.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_int* def_int(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define an integer variable with a specified key, array, shape, and info.
     * @param key The key for the variable.
     * @param arr_in Pointer to the input array.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_int* def_int(const std::string& key, kids_int* arr_in, Shape S = 1, const std::string& info = "");

    /**
     * Define an integer variable with a specified key, reference key, and info.
     * @param key The key for the variable.
     * @param key_in The key of the variable to reference.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_int* def_int(const std::string& key, const std::string& key_in, const std::string& info = "");

    /**
     * Define an integer variable with a specified key, shape, and info and update the dataset.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Reference to the updated DataSet.
     */
    DataSet& _def_int(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define a real variable with a specified key, shape, and info.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_real* def_real(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define a real variable with a specified key, array, shape, and info.
     * @param key The key for the variable.
     * @param arr_in Pointer to the input array.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_real* def_real(const std::string& key, kids_real* arr_in, Shape S = 1, const std::string& info = "");

    /**
     * Define a real variable with a specified key, reference key, and info.
     * @param key The key for the variable.
     * @param key_in The key of the variable to reference.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_real* def_real(const std::string& key, const std::string& key_in, const std::string& info = "");

    /**
     * Define a real variable with a specified key, shape, and info and update the dataset.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Reference to the updated DataSet.
     */
    DataSet& _def_real(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define a complex variable with a specified key, shape, and info.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_complex* def_complex(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define a complex variable with a specified key, array, shape, and info.
     * @param key The key for the variable.
     * @param arr_in Pointer to the input array.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_complex* def_complex(const std::string& key, kids_complex* arr_in, Shape S = 1, const std::string& info = "");

    /**
     * Define a complex variable with a specified key, reference key, and info.
     * @param key The key for the variable.
     * @param key_in The key of the variable to reference.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    kids_complex* def_complex(const std::string& key, const std::string& key_in, const std::string& info = "");

    /**
     * Define a complex variable with a specified key, shape, and info and update the dataset.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Reference to the updated DataSet.
     */
    DataSet& _def_complex(const std::string& key, Shape S = 1, const std::string& info = "");

    /**
     * Define a variable with a specified key and reference key, and info.
     * @param key The key for the variable.
     * @param key_in The key of the variable to reference.
     * @param info Additional information about the variable.
     * @return Reference to the updated DataSet.
     */
    DataSet& _def(const std::string& key, const std::string& key_in, const std::string& info = "");

    /**
     * Undefine a variable with a specified key.
     * @param key The key for the variable to undefine.
     * @return Reference to the updated DataSet.
     */
    DataSet& _undef(const std::string& key);

    /**
     * Obtain information about a variable with a specified key.
     * @param key The key for the variable.
     * @return A tuple containing the data type, pointer to data, and shape of the variable.
     */
    std::tuple<kids_dtype, void*, Shape*> obtain(const std::string& key);

    /**
     * Inquiry a specified key.
     * @param key The key for the variable.
     * @return Bool whether the key exists
     */
    bool haskey(const std::string& key);

    /**
     * Get the node corresponding to a variable with a specified key.
     * @param key The key for the variable.
     * @return Pointer to the corresponding Node.
     */
    Node* node(const std::string& key);

    /**
     * Access the DataSet corresponding to a variable with a specified key.
     * @param key The key for the variable.
     * @return Pointer to the corresponding DataSet.
     */
    DataSet* at(const std::string& key);

    /**
     * Get help information for a variable with a specified name.
     * @param name The name of the variable.
     * @return Help information as a string.
     */
    virtual std::string help(const std::string& name);

    /**
     * Get a string representation of the DataSet.
     * @return String representation of the DataSet.
     */
    virtual std::string repr();

    /**
     * Dump the DataSet to an output stream.
     * @param os The output stream to dump the DataSet (fixed in field) to.
     */
    virtual void dump_match(std::ostream& os, const std::string& prefix);

    /**
     * Dump the DataSet to an output stream.
     * @param os The output stream to dump the DataSet to.
     */
    virtual void dump(std::ostream& os);

    /**
     * Load the DataSet from an input stream.
     * @param is The input stream to load the DataSet from.
     */
    virtual void load(std::istream& is);

   private:
    /**
     * Define a variable of type T with a specified key, shape, and info.
     * @tparam T The type of variable to define.
     * @param key The key for the variable.
     * @param S The shape of the variable.
     * @param info Additional information about the variable.
     * @return Pointer to the defined variable.
     */
    template <typename T>
    T* def(const std::string& key, Shape S = 1, const std::string& info = "");

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
};

};  // namespace PROJECT_NS

#endif  // KIDS_DataSet_H

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

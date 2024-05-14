/**@file        Kernel.h
 * @brief       this file provide Kernel class
 * @details
 *  Kernal is one of the most important classes in KIDS, and is responsible
 *  for the implementment of algorithms and functionalities
 * ## Easier interfaces for getting parameters (@todo)
 *  - Dependencies: Param & DataSet. The parameters are passed by Param class.
 *    It is also free of storage of data, which is mapped to DataSet class.
 *  - It's a (hierarchical) tree structure, each kernal also consists of a
 *    list of other kernels.
 *  - The detailed algorithms can be realized in itself or its child. In
 *    priciple, we don't recomment multilevel inheritance.
 *    (please use `final` keyword).
 *  - There should be no `new` or `delete` in Kernel and its derivatives!
 *    Be clear that there are three types attributes for members in Kernel:
 *    [TODO]
 *    a) external: means only need kernel's pointer points to correct object.
 *    b) shared: means kernel's pointer points to correct object and takes
 *            ownership.
 *    c) internal: means a non-pointer data type which is not exposed to
 *              public domain.
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
 * @par [logs]:
 * <table>
 * <tr><th> Date        <th> Description
 * <tr><td> 2024-04-14  <td> Update the file.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_Kernel_H
#define KIDS_Kernel_H

#include <memory>

#include "kids/DataSet.h"
#include "kids/Param.h"
#include "kids/RuleSet.h"
#include "kids/Status.h"

namespace PROJECT_NS {

/**
 * this class provides the container and implementation of algorithms
 */
class Kernel : public std::enable_shared_from_this<Kernel> {
   public:
    /**
     * @brief Get the name of the kernel.
     *
     * Returns the concatenated name of the kernel, prefixed with "Kernel__".
     *
     * @return The name of the kernel.
     */
    inline virtual const std::string getName() { return utils::concat("Kernel__", kernel_name); }

    static int dump(Kernel* ker) { return 0; }

    static int report(Kernel* ker) { return 0; }

    static int mpi_reduce_info(Kernel* ker) { return 0; }

    /**
     * Constructor with an optional specified name.
     *
     * @param customized_name The name of the kernel (optional).
     */
    Kernel(const std::string& customized_name = "");

    /**
     * Destructor.
     */
    virtual ~Kernel();

    /**
     * Set input parameters for the kernel and its children.
     *
     * @param PM Shared pointer to the Param object.
     */
    void setInputParam(std::shared_ptr<Param>& PM);

    /**
     * Set input data set for the kernel and its children.
     *
     * @param DS Shared pointer to the DataSet object.
     */
    void setInputDataSet(std::shared_ptr<DataSet>& DS);

    /**
     * Get the parameter associated with the kernel.
     *
     * @return Shared pointer to the Param object.
     */
    std::shared_ptr<Param> getParam() const;

    /**
     * Get the data set associated with the kernel.
     *
     * @return Shared pointer to the DataSet object.
     */
    std::shared_ptr<DataSet> getDataSet() const;

    /**
     * Prepare initial conditions for the kernel and its children.
     *
     * @param stat The status object to store initialization status.
     * @return The status object after initialization.
     */
    Status& initializeKernel(Status& stat);

    /**
     * Execute the kernel's algorithm and those of its children.
     *
     * @param stat The status object to store execution status.
     * @return The status object after execution.
     */
    Status& executeKernel(Status& stat);

    /**
     * Finalize the kernel and its children, performing any necessary cleanup.
     *
     * @param stat The status object to store finalization status.
     * @return The status object after finalization.
     */
    Status& finalizeKernel(Status& stat);

    /**
     * Get the type of the kernel.
     *
     * This function returns the type of the kernel.
     *
     * @note This function should be implemented in derived classes to return the specific type of the kernel.
     *
     * @return The type of the kernel.
     */
    virtual int getType() const;

    /**
     * Get the ID of the kernel.
     *
     * @return The ID of the kernel.
     */
    int getID() const;

    /**
     * Overloaded equality operator to compare two Kernel objects by their IDs.
     *
     * @param ker The Kernel object to compare with.
     * @return True if the IDs of the two kernels are equal, otherwise false.
     */
    bool operator==(const Kernel& ker);

    /**
     * Append a kernel as the last child of the current tree node.
     *
     * @param ker The kernel to append as the last child.
     * @return Reference to the modified kernel.
     */
    Kernel& appendChild(std::shared_ptr<Kernel> ker);

    /**
     * Insert a kernel at specified indexes in the tree.
     *
     * @param indexes Indexes indicating the position to insert the kernel.
     * @param ker The kernel to insert.
     * @return Reference to the modified kernel.
     */
    Kernel& insertAt(std::vector<std::size_t> indexes, std::shared_ptr<Kernel> ker);

    /**
     * Remove kernels at specified indexes from the tree.
     *
     * @param indexes Indexes indicating the kernels to remove.
     * @return Reference to the modified kernel.
     */
    Kernel& removeAt(std::vector<std::size_t> indexes);

    /**
     * Update the kernel at specified indexes in the tree.
     *
     * @param indexes Indexes indicating the kernel to update.
     * @param ker The kernel to update with.
     * @return Reference to the modified kernel.
     */
    Kernel& updateAt(std::vector<std::size_t> indexes, std::shared_ptr<Kernel> ker);

    /**
     * Retrieve the last parent kernel along with the order of its child kernels, if available.
     *
     * @return A tuple containing a pointer to the last parent kernel and the order of its child kernels.
     */
    std::tuple<Kernel*, std::size_t> getLastParentKernelAndChildOrder();

    /**
     * Get RuleSet associated with the Kernel
     *
     * @return     The RuleSet (name, dataPointer, datatype, size, stride).
     */
    std::shared_ptr<RuleSet> getRuleSet();

    /**
     * Serialize a Kernel object into a string representation.
     *
     * @param ker The Kernel object to serialize.
     * @return A string representing the serialized Kernel object.
     */
    static const std::string serializeKernel(const Kernel& ker);

    /**
     * Deserialize a string representation into a Kernel object.
     *
     * @param str The string containing the serialized Kernel object.
     * @return A shared pointer to the deserialized Kernel object.
     */
    static std::shared_ptr<Kernel> deserializeKernel(const std::string& str);

    /**
     * Generate a formatted string containing information about the kernel.
     *
     * This function generates a formatted string containing information about the kernel,
     * including its total time, current layer, total depth, and total alignment size.
     *
     * @param total_time The total time taken by the kernel (default: -1.0f).
     * @param current_layer The current layer of the kernel (default: 0).
     * @param total_depth The total depth of the kernel (default: 0).
     * @param total_align_size The total alignment size of the kernel (default: 0).
     * @return A formatted string containing kernel information.
     */
    const std::string generateInformationString(double total_time       = -1.0f,  //
                                                int    current_layer    = 0,      //
                                                int    total_depth      = 0,      //
                                                int    total_align_size = 0);


   protected:
    /**
     * @defgroup KernelMetadata Metadata for the Kernel class
     * @{
     */

    bool        is_timing      = false;  ///< Flag indicating whether timing is enabled for this kernel.
    bool        has_parent     = false;  ///< Flag indicating whether the kernel has a parent.
    int         count_calc     = 0;      ///< Counter for the number of calculations performed by this kernel.
    int         count_exec     = 0;      ///< Counter for the number of executions performed by this kernel.
    int         kernel_id      = 0;      ///< ID of the kernel.
    int         kernel_type    = 0;      ///< Type of the kernel.
    double      exec_time      = 0.0f;   ///< Total execution time of the kernel.
    int         depth          = 0;      ///< Depth of the kernel in the tree structure.
    int         max_align_size = 0;      ///< Maximum alignment size used by this kernel.
    std::string kernel_name;             ///< Name of the kernel.

    /**
     * @}
     */

    /**
     * @brief Shared pointer to the Param object associated with this kernel.
     */
    std::shared_ptr<Param> _param;

    /**
     * @brief Shared pointer to the DataSet object associated with this kernel.
     */
    std::shared_ptr<DataSet> _dataset;

    /**
     * @brief Pointer to the parent kernel.
     */
    Kernel* _parent_kernel;

    /**
     * @brief Order of this kernel in its parent's children.
     */
    std::size_t _order_in_parent;

    /**
     * @brief Vector containing shared pointers to the child kernels of this kernel.
     */
    std::vector<std::shared_ptr<Kernel>> _child_kernels;

    /**
     * @brief Vector containing shared pointers to all descendant kernels of this kernel.
     */
    std::vector<std::shared_ptr<Kernel>> _all_kernels;

    /**
     * @brief Recorded Rules associated with the Kernel.
     */
    std::shared_ptr<RuleSet> _ruleset;

    /**
     * @brief Virtual function to set input parameters for the kernel implementation.
     *
     * @param PM Shared pointer to the Param object containing input parameters.
     */
    virtual void setInputParam_impl(std::shared_ptr<Param>& PM);

    /**
     * @brief Virtual function to set input data set for the kernel implementation.
     *
     * @param DS Shared pointer to the DataSet object containing input data.
     */
    virtual void setInputDataSet_impl(std::shared_ptr<DataSet>& DS);

    /**
     * @brief Virtual function to initialize the kernel implementation.
     *
     * @param stat Status object to store initialization status.
     * @return Status object after initialization.
     */
    virtual Status& initializeKernel_impl(Status& stat);

    /**
     * @brief Virtual function to execute the kernel implementation.
     *
     * @param stat Status object to store execution status.
     * @return Status object after execution.
     */
    virtual Status& executeKernel_impl(Status& stat);

    /**
     * @brief Virtual function to finalize the kernel implementation.
     *
     * @param stat Status object to store finalization status.
     * @return Status object after finalization.
     */
    virtual Status& finalizeKernel_impl(Status& stat);

   private:
    /**
     * @brief Connect related kernels to this kernel.
     *
     * This function connects related kernels to the current kernel.
     * Related kernels are kernels that are associated with or connected to this kernel in some way.
     *
     * @note The exact definition of "related kernels" may vary depending on the context and application.
     */
    void connectRelatedKernels(std::shared_ptr<Kernel>& ker);

    /**
     * @brief Get the dictionary of kernels (mapping from names to kernel pointers).
     *
     * This function returns a reference to the static map containing the dictionary of kernels.
     * The dictionary maps kernel names to their corresponding kernel pointers.
     *
     * @return A reference to the dictionary of kernels.
     */
    static std::map<std::string, Kernel*>& getDictOfKernels();

    /**
     * @brief Get the vector of all kernel pointers.
     *
     * This function returns a reference to the static vector containing pointers to all kernels.
     *
     * @return A reference to the vector of all kernel pointers.
     */
    static std::vector<Kernel*>& getKernels();
};

};  // namespace PROJECT_NS

#endif  // KIDS_Kernel_H

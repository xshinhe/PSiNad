/**@file        Einsum.h
 * @brief       this file provides einsum operation
 * @details
 * ## About the rules in einsum operations provided by "einsum.h"
 *
 * In this documentation, we define key terms related to the einsum operation:
 * - einsum label: Refers to the labels assigned to the indices of the operands
 *   involved in the einsum operation. For example, 'a', 'b', 'c' could be operand
 *   labels.
 * - einsum shape: Describes the shape or dimensions of a tensor corresponding to
 *   each operand. For instance, "abc" could represent the shape of a tensor.
 * - einsum expression: Specifies the einsum operation to be performed, indicating
 *   the contraction pattern of indices between operands. For example, "ab,bc->ac"
 *   defines how the indices from two operands are contracted to produce the final
 *   output.
 *
 * einsum label:
 * - fixed label: expressed as `[name]`. It cannot appear in user's deduction!
 * - inner label: appeared in user's deduction or only appears once in summation.
 * - outer label: not in user's deduction or appears multiple times in summation.
 *
 * rules for einsum expression:
 * - user's deduction: `->` indicates for einsum result. Then all (non-fixed)
 *   labels appear behind `->` are treated as outer labels, while those not appear
 *   behind `->` are used as inner labels.
 * - automatic deduction: if `->` is not provided in einsum expression. The inner/
 *   outer labels are seperated by the times they appear in summation.
 *
 * ## Declaration
 *
 * This file provides a einsum operation with an easy realization, the performance
 * is not optimized so well. Based on some test, there are about 20x slower than
 * in optimized loop in pure C.
 *
 * ## benckmark cases
 *
 * the following benchmark test is passed:
 * ```cpp
 * std::size_t L = 2 * 2 * 2 * 2 * 3 * 3 * 3 * 3;
 * std::vector<int> A(L, 0);
 * for (int z = 0; z < L; ++z) {
 *     A[z] = z % 3 - (z % 5) * (z % 5) + (z % 7);
 * }
 * std::vector<int> res(L, 0);
 * einsum("i", {A.data()}, {{L}}, res.data(), {L});
 * ARRAY_SHOW(res.data(), 1, 1);
 * einsum("i->", {A.data()}, {{L}}, res.data(), {1});
 * ARRAY_SHOW(res.data(), 1, 1);
 * einsum("ikkkji->j", {A.data()}, {{2, 3, 3, 3, 12, 2}}, res.data(), {12});
 * ARRAY_SHOW(res.data(), 1, 12);
 * einsum("ikkkji->ik", {A.data()}, {{2, 3, 3, 3, 12, 2}}, res.data(), {2, 3});
 * ARRAY_SHOW(res.data(), 2, 3);
 * einsum("i,i", {A.data(), A.data()}, {{L}, {L}}, res.data(), {1});
 * ARRAY_SHOW(res.data(), 1, 1);
 * einsum("ik,ik", {A.data(), A.data()}, {{16, 81}, {16, 81}}, res.data(), {1});
 * ARRAY_SHOW(res.data(), 1, 1);
 * einsum("ik,ki", {A.data(), A.data()}, {{16, 81}, {81, 16}}, res.data(), {1});
 * ARRAY_SHOW(res.data(), 1, 1);
 * einsum("ik,kj->ij", {A.data(), A.data()}, {{4, 324}, {324, 4}},
 *                     res.data(), {4, 4});
 * ARRAY_SHOW(res.data(), 4, 4);
 * einsum("ik,kj,ljjlll->il", {A.data(), A.data(), A.data()},
 *        {{4, 324}, {324, 4}, {3, 4, 4, 3, 3, 3}},
 *        res.data(), {4, 3});
 * ARRAY_SHOW(res.data(), 4, 3);
 * ```
 * which is coincided with the results by numpy's einsum() function:
 * ```py
 * import numpy as np
 * L = 2*2*2*2*3*3*3*3
 * z = np.arange(L)
 * A = (z%3) - (z%5)**2 + (z%7)
 * print(np.einsum('i', A))
 * print(np.einsum('i->', A))
 * print(np.einsum('ikkkji->j', A.reshape((2,3,3,3,12,2))))
 * print(np.einsum('ikkkji->ik', A.reshape((2,3,3,3,12,2))))
 * print(np.einsum('i,i', A, A))
 * print(np.einsum('ik,ik', A.reshape((16,81)), A.reshape((16,81))))
 * print(np.einsum('ik,ki', A.reshape((16,81)), A.reshape((81,16))))
 * print(np.einsum('ik,kj->ij', A.reshape((4,324)), A.reshape((324,4))))
 * print(np.einsum('ik,kj,ljjlll->il', A.reshape((4,324)),
 *                  A.reshape((324,4)), A.reshape((3,4,4,3,3,3))))
 * ```
 * the results are:
 * --------------------------------
 * ```text
 * [ 0  1  0 ... -4 -9  2]
 * -2589
 * [-25  -9  -9  -8 -11 -22  -8 -15  -5 -10 -21  -5]
 * [[-18 -28 -33]
 *  [-27 -21 -21]]
 * 56293
 * 56293
 * 52891
 * [[-2724  3071  8230  7447]
 *  [-8259 -3074  3047  8535]
 *  [ 6290 -8694 -3109  3442]
 *  [ 8340  6031 -8557 -2810]]
 * [[ 83744  54125  79687]
 *  [ 57250 -22579  -1748]
 *  [ -9172 -21858 -39479]
 *  [-39026  55563 -10344]]
 * ```
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
 * <tr><td> 2024-03-29  <td> initial version. The shape_inputs and shape_output are
 *                          consistent with that interface Basic_Shape in the file
 *                          Shape.h.
 * </table>
 *
 **********************************************************************************
 */

#ifndef KIDS_EINSUM_H
#define KIDS_EINSUM_H

#include <cstring>
#include <string>
#include <vector>

namespace PROJECT_NS {

//*********************************************************************************
/**
 * EinsumIdx is a struct store information of index used in einsum operation
 */
struct EinsumIdx {
    EinsumIdx(char _label, std::size_t _cnt, std::size_t _dim, std::size_t _val)
        : label{_label}, cnt{_cnt}, dim{_dim}, val{_val} {};

    char        label;    ///< unique identifer for EinsumIdx
    std::size_t cnt = 0;  ///< indicate the type (0: fixed; 1: outer; >1: inner)
    std::size_t dim = 0;  ///< bound of the value of index
    std::size_t val = 0;  ///< the value of the index
};

/**
 * DimenHelper is a struct control dimensional utils on the orginal/einsum index
 * for a given tensor.
 */
struct DimenHelper {
   public:
    std::size_t esshape_rank;  ///< the rank of the tensor
    std::size_t total_esidx;   ///< size if the EinsumIdx System
    std::size_t total_size;
    std::size_t fixed_init_idx;

    std::vector<std::size_t> dims;      ///< leading dimensions of the tensor
    std::vector<std::size_t> ldims;     ///< leading dimensions of the tensor
    std::vector<std::size_t> es_ldims;  ///< leading dimensions of the tensor represented in einsum indexes
    std::vector<std::size_t> mapldims;  ///< utils for sum of several leading dimensions as the shift step

    DimenHelper(){};

    /**
     * @param[in]       esshape      literals describing the rule of a tensor like
     *                                  i.e. "abc" represents a rank-3 tensor
     * @param[in]       idx_vec         EinsumIdx System (std::vector<EinsumIdx>&)
     */
    DimenHelper(const std::string& esshape, std::vector<EinsumIdx>& idx_vec);
};


class EinsumHelper {
   public:
    std::size_t total_esidx;   ///< total number of EinsumIdx in EinsumIdx System
    std::size_t total_tensor;  ///< total number of tensor in einsum rule

    std::vector<EinsumIdx>   einsum_idxs;  ///< the EinsumIdx System
    std::vector<std::size_t> einsum_dims;  ///< each dimension of EinsumIdx System

    std::vector<std::string> fixed_label_names;  ///< store for fixed labels

    std::vector<std::string> esshape_inputs;       ///< store einsum's strings of input tensors
    std::string              esshape_output = "";  ///< store/deduct einsum's for the ouput tensor

    std::vector<DimenHelper> dh_inputs;  ///< DimenHelper for input tensors
    DimenHelper              dh_output;  ///< DimenHelper for ouput tensor

    std::vector<std::size_t> einsum_iposes;     ///< idx placeholder for EinsumIdx System
    std::vector<std::size_t> ipos_inputs_init;  ///< idx placeholder for inital input tensors
    std::vector<std::size_t> ipos_inputs;       ///< idx placeholder for input tensors

    int count1     = 0;
    int count2     = 0;
    int count3     = 0;
    int total_loop = 0;

    /**
     * @param[in]       einsum_expression      expression for einsum rule
     * @param[in]       shape_inputs           input shapes as a vector
     * @param[in]       shape_output           output shapes
     */
    EinsumHelper(const std::string&                    einsum_expression,  //
                 std::vector<std::vector<std::size_t>> shape_inputs,       //
                 std::vector<std::size_t>              shape_output = {}   //
    );
};

/**
 * @tparam          T               data type
 * @param[in]       EH              EinsumHelper object
 * @param[in]       data_inputs     vector of pointers of data of input tensors
 * @param[inout]    data_output     pointer stored data of output tensor
 */
template <typename T>
void einsum(EinsumHelper&          EH,           //
            const std::vector<T*>& data_inputs,  //
            T*                     data_output   //
) {
    auto& einsum_dims      = EH.einsum_dims;
    auto& einsum_iposes    = EH.einsum_iposes;
    auto& ipos_inputs      = EH.ipos_inputs;
    auto& ipos_inputs_init = EH.ipos_inputs_init;
    // ipos_output
    auto& dh_inputs          = EH.dh_inputs;
    auto& dh_output_mapldims = EH.dh_output.mapldims;

    std::size_t total_loop   = EH.total_loop;
    std::size_t total_tensor = EH.total_tensor;
    std::size_t total_esidx  = EH.total_esidx;
    std::size_t imax         = EH.count3 - 1;
    std::size_t imin         = EH.count1;

    memset(einsum_iposes.data(), 0, total_esidx * sizeof(std::size_t));
    memcpy(ipos_inputs.data(), ipos_inputs_init.data(), total_tensor * sizeof(std::size_t));
    bool reset_zero = true;
    data_output[0]  = T(0);
    for (std::size_t iloop = 0, ipos_output = 0; iloop < total_loop; ++iloop) {
        if (reset_zero) data_output[ipos_output] = T(0);

        T term = T(1);
        for (int iten = 0; iten < total_tensor; ++iten) { term *= data_inputs[iten][ipos_inputs[iten]]; }
        data_output[ipos_output] += term;

        std::size_t i = imax;
        while (++einsum_iposes[i] == einsum_dims[i] && i > imin) { einsum_iposes[i--] = 0; }
        reset_zero = (i < EH.count2);

        for (int iten = 0; iten < total_tensor; ++iten)  //
            ipos_inputs[iten] += dh_inputs[iten].mapldims[i];
        ipos_output += dh_output_mapldims[i];
    }
}


/**
 * @tparam          T                   data type
 * @param[in]       einsum_expression   expression for einsum rule
 * @param[in]       data_inputs         vector of pointers of data of input tensors
 * @param[in]       shapes_inputs       vector of shapes of input tensors
 * @param[inout]    data_output         pointer stored data of output tensor
 * @param[in]       shape_output        shape of the output tensor
 */
template <typename T>
void einsum(const std::string&                           einsum_expression,  //
            std::vector<T*>                              data_inputs,        //
            const std::vector<std::vector<std::size_t>>& shape_inputs,       //
            T*                                           data_output,        //
            const std::vector<std::size_t>&              shape_output = {}   //
) {
    EinsumHelper EH(einsum_expression, shape_inputs, shape_output);
    einsum(EH, data_inputs, data_output);
}

};  // namespace PROJECT_NS

#endif  // KIDS_EINSUM_H

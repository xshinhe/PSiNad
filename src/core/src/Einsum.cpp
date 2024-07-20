#include "kids/Einsum.h"

#include <algorithm>
#include <iostream>
#include <sstream>

namespace PROJECT_NS {

DimenHelper::DimenHelper(const std::string& esshape, std::vector<EinsumIdx>& idx_vec)
    : esshape_rank{esshape.size()}, total_esidx{idx_vec.size()} {
    dims.resize(esshape_rank);
    ldims.resize(esshape_rank);
    es_ldims.resize(total_esidx);
    mapldims.resize(total_esidx);

    // calculate the normal leading dimensions of the tensor
    for (int k = esshape_rank - 1, lastsize = 1, lastldim = 1; k >= 0; --k) {
        ldims[k] = lastsize * lastldim;
        int q    = -1;
        while (idx_vec[++q].label != esshape[k]) {};
        dims[k]  = idx_vec[q].dim;
        lastsize = dims[k];
        lastldim = ldims[k];
    }

    total_size = 1;
    for (auto& i : dims) total_size *= i;

    // calculate the leading dimensions of the tensor represented in einsum indexes
    for (int i = 0; i < total_esidx; ++i) {
        char c      = idx_vec[i].label;
        es_ldims[i] = 0;
        for (int k = esshape_rank - 1; k >= 0; --k) {
            if (c == esshape[k]) es_ldims[i] += ldims[k];
        }
    }

    // calculate sum of several leading dimensions as the shift step
    for (int i = 0; i < total_esidx; ++i) {  //
        mapldims[i] = es_ldims[i];
        for (int k = i + 1; k < total_esidx; ++k) {  //
            mapldims[i] -= (idx_vec[k].dim - 1) * es_ldims[k];
        }
    }
}

EinsumHelper::EinsumHelper(const std::string&                    einsum_expression,  //
                           std::vector<std::vector<std::size_t>> shape_inputs,       //
                           std::vector<std::size_t>              shape_output        //
) {
    std::stringstream ss{einsum_expression};
    std::string       esshape        = "";
    int               ishape         = 0;
    bool              auto_deduction = true;
    for (char c; ss >> c;) {
        switch (c) {
            case ',': {
                if (esshape.size() != shape_inputs[ishape].size()) {
                    std::cerr << "esshape.size() = " << esshape.size() << "\n";
                    std::cerr << "shape_inputs[ishape].size() = " << shape_inputs[ishape].size() << "\n";
                    throw std::runtime_error("mismatch einsum rule with shape!");
                }
                esshape_inputs.push_back(esshape);
                esshape = "";
                ishape++;
                break;
            }
            case '[': {
                std::string label_name = "";
                while (ss >> c) {
                    if (c == ']') break;
                    label_name += c;
                }
                auto it    = std::find(fixed_label_names.begin(), fixed_label_names.end(), label_name);
                auto found = (it != fixed_label_names.end());
                int  ipos  = found ? int(it - fixed_label_names.begin()) : fixed_label_names.size();
                c          = (char) ((int) '0' + ipos);

                if (!found) {
                    einsum_idxs.push_back(EinsumIdx(c, 0, shape_inputs[ishape][esshape.size()], 0));
                    fixed_label_names.push_back(label_name);
                } else {
                    auto it2 = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),
                                            [c](EinsumIdx idx) { return c == idx.label; });
                    if (it2->dim != shape_inputs[ishape][esshape.size()]) {
                        // std::cout << c << shape_inputs[ishape][esshape.size()] << "\n";
                        throw std::runtime_error("bad einsum shape!");
                    }
                }
                esshape += c;

                if (fixed_label_names.size() > 10) throw std::runtime_error("too many fixed einsum idx!");
                break;
            }
            case ' ':
            case '-':
                break;
            case '>': {
                auto_deduction = false;  // then by user's deduction
                esshape_output = "";
                while (ss >> c) {
                    if ((int) c < (int) 'a' || (int) c > (int) 'z') {
                        throw std::runtime_error("only allowed [a-z] for normal einsum label");
                    }
                    auto it = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),  //
                                           [c](EinsumIdx idx) { return c == idx.label; });
                    if (it != einsum_idxs.end()) {
                        if (shape_output.size() > 0 && it->dim != shape_output[esshape_output.size()]) {
                            throw std::runtime_error("bad einsum shape!");
                        }
                        it->cnt = 1;
                    } else {
                        std::cerr << esshape_output << "\n";
                        throw std::runtime_error("bad einsum einsum_expression!");
                    }
                    esshape_output += c;
                }
                if (esshape_output.size() != shape_output.size() && shape_output.size() != 0) {
                    // if shape_output.size() == 0, we don't check
                    throw std::runtime_error("mismatch einsum rule with shape!");
                }
                break;
            }
            default: {
                if ((int) c < (int) 'a' || (int) c > (int) 'z') {
                    throw std::runtime_error("only allowed [a-z] for normal einsum label");
                }
                auto it = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),  //
                                       [c](EinsumIdx idx) { return c == idx.label; });
                if (it != einsum_idxs.end()) {
                    if (it->dim == shape_inputs[ishape][esshape.size()]) {
                        it->cnt++;  // update as inner label
                    } else {
                        // std::cout << c << shape_inputs[ishape][esshape.size()] << "\n";
                        throw std::runtime_error("bad einsum shape!");
                    }
                } else {
                    einsum_idxs.push_back(EinsumIdx(c, 1, shape_inputs[ishape][esshape.size()], 0));
                }
                esshape += c;
                break;
            }
        }
    }
    esshape_inputs.push_back(esshape);

    if (auto_deduction) {
        esshape_output = "";
        for (auto& idx : einsum_idxs) {
            if (idx.cnt == 1) esshape_output += idx.label;
        }
    } else {
        for (auto& idx : einsum_idxs) {
            if (idx.cnt == 1) idx.cnt = 2;  // revise to inner label
            for (auto& label : esshape_output) {
                if (idx.label == label) idx.cnt = 1;  // revise to outer label
            }
        }
    }
    if (esshape_output == "") {  // allow return a scalar
        einsum_idxs.push_back(EinsumIdx('*', 0, 1, 0));
        esshape_output = "*";
    }
    std::sort(einsum_idxs.begin(), einsum_idxs.end(),
              [](EinsumIdx idx1, EinsumIdx idx2) { return idx1.cnt < idx2.cnt; });

    count1     = 0;
    count2     = 0;
    count3     = 0;
    total_loop = 1;
    for (auto& idx : einsum_idxs) {
        if (idx.cnt <= 0) count1++;
        if (idx.cnt <= 1) count2++;
        if (idx.cnt > 0) total_loop *= idx.dim;
        count3++;
    }

    // for (auto& idx : einsum_idxs) {
    //     std::cout << idx.label << ", " << idx.cnt << "," << idx.dim << ", " << idx.val << "\n";
    // }
    // for (auto& esshape : esshape_inputs) std::cout << esshape << "\n";
    // std::cout << "->" << esshape_output << "\n";
    // std::cout << "count1 : " << count1 << "\n";
    // std::cout << "count2 : " << count2 << "\n";
    // std::cout << "count3 : " << count3 << "\n";

    for (auto& esshape : esshape_inputs) { dh_inputs.push_back(DimenHelper(esshape, einsum_idxs)); }
    dh_output = DimenHelper(esshape_output, einsum_idxs);

    for (auto& idx : einsum_idxs) { einsum_dims.push_back(idx.dim); }

    total_esidx  = einsum_idxs.size();
    total_tensor = esshape_inputs.size();

    einsum_iposes.resize(total_esidx);
    ipos_inputs.resize(total_tensor);
}

};  // namespace PROJECT_NS
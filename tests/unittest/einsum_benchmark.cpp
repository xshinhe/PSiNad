#include <algorithm>
#include <chrono>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

/**
 * control the output printing format
 */
constexpr inline int FMT_WIDTH(int X) { return X + 7; }
#define FMT(X)                                                            \
    " " << std::setiosflags(std::ios::scientific) /*scientific notation*/ \
        << std::setprecision(X)                   /*precision*/           \
        << std::right                             /*alignment*/           \
        << std::setw(FMT_WIDTH(X))                /*width of text*/

#define ARRAY_SHOW(_A, _n1, _n2)                                                     \
    ({                                                                               \
        std::cout << "Show Array <" << #_A << ">\n";                                 \
        int _idxA = 0;                                                               \
        for (int _i = 0; _i < (_n1); ++_i) {                                         \
            for (int _j = 0; _j < (_n2); ++_j) std::cout << FMT(0) << (_A)[_idxA++]; \
            std::cout << std::endl;                                                  \
        }                                                                            \
    })

struct EinsumIdx {
    char label;
    std::size_t cnt = 0;
    std::size_t dim = 0;
    std::size_t val = 0;
};

struct DimHelper {
   public:
    std::size_t tensor_rank;
    std::size_t total_esidx;

    std::vector<std::size_t> ldims;
    std::vector<std::size_t> es_ldims;
    std::vector<std::size_t> mapldims;

    DimHelper(){};

    DimHelper(const std::string& tensor_str, std::vector<EinsumIdx>& idx_vec)
        : tensor_rank{tensor_str.size()}, total_esidx{idx_vec.size()} {
        ldims.resize(tensor_rank);
        es_ldims.resize(total_esidx);
        mapldims.resize(total_esidx);

        // get normal ldims ( size = tensor_rank )
        for (int k = tensor_rank - 1, lastsize = 1, lastldim = 1; k >= 0; --k) {
            ldims[k] = lastsize * lastldim;
            int q    = -1;
            while (idx_vec[++q].label != tensor_str[k]) {};
            lastsize = idx_vec[q].dim;
            lastldim = ldims[k];
        }

        // get einsum ldims ( size = total_esidx )
        for (int i = 0; i < total_esidx; ++i) {
            char c      = idx_vec[i].label;
            es_ldims[i] = 0;
            for (int k = tensor_rank - 1; k >= 0; --k) {
                if (c == tensor_str[k]) es_ldims[i] += ldims[k];
            }
        }

        // get mapldims
        for (int i = 0; i < total_esidx; ++i) {  //
            mapldims[i] = es_ldims[i];
            for (int k = i + 1; k < total_esidx; ++k) {  //
                mapldims[i] -= (idx_vec[k].dim - 1) * es_ldims[k];
            }
        }
    }
};


class EinsumHelper {
   public:
    std::size_t total_esidx;
    std::size_t total_tensor;

    std::vector<EinsumIdx> einsum_idxs;
    std::vector<std::size_t> einsum_dims;

    std::vector<std::string> fixed_exprs;

    std::vector<std::string> tensor_str_inputs;
    std::string tensor_str_output = "";

    std::vector<DimHelper> dh_inputs;
    DimHelper dh_output;

    std::vector<std::size_t> einsum_iposes;
    std::vector<std::size_t> ipos_inputs;

    int count1     = 0;
    int count2     = 0;
    int count3     = 0;
    int total_loop = 0;

    EinsumHelper(const std::string& einsum_rule,                      //
                 std::vector<std::vector<std::size_t>> shape_inputs,  //
                 std::vector<std::size_t> shape_output = {}           //
    ) {
        std::stringstream ss{einsum_rule};
        std::string tensor_str = "";
        int ishape             = 0;
        for (char c; ss >> c;) {
            switch (c) {
                case ',': {
                    tensor_str_inputs.push_back(tensor_str);
                    tensor_str = "";
                    ishape++;
                    break;
                }
                case '[': {
                    std::string fixedname = "";
                    while (ss >> c) {
                        if (c == ']') break;
                        fixedname += c;
                    }
                    auto it    = std::find(fixed_exprs.begin(), fixed_exprs.end(), fixedname);
                    auto found = (it != fixed_exprs.end());
                    int ipos   = found ? int(it - fixed_exprs.begin()) : fixed_exprs.size();
                    c          = (char) ((int) '0' + ipos);

                    if (!found) {
                        einsum_idxs.push_back(EinsumIdx{.label = c,
                                                        .cnt   = 0,  //
                                                        .dim   = shape_inputs[ishape][tensor_str.size()],
                                                        .val   = 0});
                        fixed_exprs.push_back(fixedname);
                    } else {
                        auto it2 = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),
                                                [c](EinsumIdx idx) { return c == idx.label; });
                        if (it2->dim != shape_inputs[ishape][tensor_str.size()]) {
                            std::cout << c << shape_inputs[ishape][tensor_str.size()] << "\n";
                            throw std::runtime_error("bad einsum shape!");
                        }
                    }
                    tensor_str += c;

                    if (fixed_exprs.size() > 10) throw std::runtime_error("too many fixed einsum idx!");
                    break;
                }
                case ' ':
                case '-':
                    break;
                case '>': {
                    tensor_str_output = "";
                    while (ss >> c) {
                        auto it = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),  //
                                               [c](EinsumIdx idx) { return c == idx.label; });
                        if (it != einsum_idxs.end()) {
                            if (shape_output.size() > 0 && it->dim != shape_output[tensor_str_output.size()]) {
                                throw std::runtime_error("bad einsum shape!");
                            }
                            it->cnt = 1;
                        } else {
                            throw std::runtime_error("bad einsum einsum_rule!");
                        }
                        tensor_str_output += c;
                    }
                    break;
                }
                default: {
                    auto it = std::find_if(einsum_idxs.begin(), einsum_idxs.end(),  //
                                           [c](EinsumIdx idx) { return c == idx.label; });
                    if (it != einsum_idxs.end()) {
                        if (it->dim == shape_inputs[ishape][tensor_str.size()]) {
                            it->cnt++;
                        } else {
                            std::cout << c << shape_inputs[ishape][tensor_str.size()] << "\n";
                            throw std::runtime_error("bad einsum shape!");
                        }
                    } else {
                        einsum_idxs.push_back(
                            EinsumIdx{.label = c, .cnt = 1, .dim = shape_inputs[ishape][tensor_str.size()], .val = 0});
                    }
                    tensor_str += c;
                    break;
                }
            }
        }
        tensor_str_inputs.push_back(tensor_str);

        if (std::all_of(einsum_idxs.begin(), einsum_idxs.end(), [](EinsumIdx idx) { return idx.cnt > 1; })) {
            einsum_idxs.push_back(EinsumIdx{.label = '0',
                                            .cnt   = 0,  //
                                            .dim   = 1,
                                            .val   = 0});
            tensor_str_output = "0";
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
        // for (auto& s : tensor_str_inputs) std::cout << s << "\n";
        // std::cout << "fixed count : " << count1 << "\n";
        // std::cout << "fixed + outer count : " << count2 << "\n";
        // std::cout << "all count : " << count3 << "\n";

        for (auto& tensor_str : tensor_str_inputs) { dh_inputs.push_back(DimHelper(tensor_str, einsum_idxs)); }
        dh_output = DimHelper(tensor_str_output, einsum_idxs);

        for (auto& idx : einsum_idxs) { einsum_dims.push_back(idx.dim); }

        total_esidx  = einsum_idxs.size();
        total_tensor = tensor_str_inputs.size();

        einsum_iposes.resize(total_esidx);
        ipos_inputs.resize(total_tensor);
    }
};

template <typename T>
void einsum(EinsumHelper& EH,                    //
            const std::vector<T*>& data_inputs,  //
            T* data_output                       //
) {
    auto& einsum_dims   = EH.einsum_dims;
    auto& einsum_iposes = EH.einsum_iposes;
    auto& ipos_inputs   = EH.ipos_inputs;
    // ipos_output
    auto& dh_inputs          = EH.dh_inputs;
    auto& dh_output_mapldims = EH.dh_output.mapldims;

    std::size_t total_loop   = EH.total_loop;
    std::size_t total_tensor = EH.total_tensor;
    std::size_t total_esidx  = EH.total_esidx;
    std::size_t imax         = EH.count3 - 1;
    std::size_t imin         = EH.count1;

    memset(einsum_iposes.data(), 0, total_esidx * sizeof(std::size_t));
    memset(ipos_inputs.data(), 0, total_tensor * sizeof(std::size_t));
    data_output[0] = T(0);
    for (std::size_t iloop = 0, ipos_output = 0; iloop < total_loop; ++iloop) {
        T term = T(1);
        for (int iten = 0; iten < total_tensor; ++iten) { term *= data_inputs[iten][ipos_inputs[iten]]; }
        data_output[ipos_output] += term;

        std::size_t i = imax;
        while (++einsum_iposes[i] == einsum_dims[i] && i > imin) { einsum_iposes[i--] = 0; }

        for (int iten = 0; iten < total_tensor; ++iten)  //
            ipos_inputs[iten] += dh_inputs[iten].mapldims[i];

        ipos_output += dh_output_mapldims[i];
        if (i < EH.count2) data_output[ipos_output] = T(0);
    }
}

template <typename T>
void einsum(const std::string& einsum_rule,                      //
            const std::vector<T*>& data_inputs,                  //
            std::vector<std::vector<std::size_t>> shape_inputs,  //
            const T* data_output,                                //
            std::vector<std::size_t> shape_output = {}           //
) {
    EinsumHelper EH(einsum_rule, shape_inputs, shape_output);
    einsum(EH, data_inputs, data_output);
}

struct ESIDX {
   public:
    char c;  // flag
    int d;   // dimension
};

struct Dhelper {
   public:
    std::vector<int> ldims;
    std::vector<int> eldims;
    std::vector<int> sm1xeldims;
    std::vector<int> mapldims;
};

Dhelper get_ldims(const std::string& flags, std::vector<ESIDX>& idx_vec) {
    Dhelper dh;
    dh.ldims.resize(flags.size());
    dh.eldims.resize(idx_vec.size());
    dh.mapldims.resize(idx_vec.size());

    // get normal ldims:
    for (int k = flags.size() - 1, lastsize = 1, lastldim = 1; k >= 0; --k) {
        dh.ldims[k] = lastsize * lastldim;
        int q       = -1;
        while (idx_vec[++q].c != flags[k]) {};
        lastsize = idx_vec[q].d;
        lastldim = dh.ldims[k];
    }

    // get extend ldims;
    for (int i = 0; i < idx_vec.size(); ++i) {
        char c = idx_vec[i].c;

        dh.eldims[i] = 0;
        for (int k = flags.size() - 1, lastsize = 1, lastldim = 1; k >= 0; --k) {
            if (c == flags[k]) dh.eldims[i] += dh.ldims[k];
        }
    }

    // get mapldims
    for (int i = 0; i < idx_vec.size(); ++i) {  //
        dh.mapldims[i] = dh.eldims[i];
        for (int k = i + 1; k < idx_vec.size(); ++k) {  //
            dh.mapldims[i] -= (idx_vec[k].d - 1) * dh.eldims[k];
        }
    }

    //
    return dh;
}

std::string gen_einsum_redstr(std::vector<std::string> str_vec) { return "bd"; }

int main() {
    int L1 = 100;
    int L2 = 99;
    int A[L1 * L2], B[L2 * L1], AB[L1 * L1], BA[L2 * L2];

    std::vector<ESIDX> idx_vec = {
        // {'a', 2},
        {'b', L2},
        {'c', L1},
        // {'d', 5},
    };

    for (int i = 0; i < L1 * L2; ++i) {
        A[i] = i + i % 3 - (i * i) % 9;
        B[i] = 20 - i % 5 - i / 2;
    }

    int Nloop = 10000;

    int res0   = 0;
    auto begin = std::chrono::steady_clock::now();
    {
        for (int loop = 0; loop < Nloop; ++loop) {
            res0 = 0;
            for (int i = 0, idxA = 0; i < L1; ++i) {
                for (int k = 0, idxB = i; k < L2; ++k, ++idxA, idxB += L1) { res0 += A[idxA] * B[idxB]; }
            }
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "timing = " << static_cast<std::chrono::duration<double>>(end - begin).count() << "\n";
    std::cout << res0 << "?\n";

    std::string Aflag = "cb";
    std::string Bflag = "bc";
    std::string Cflag = gen_einsum_redstr({Aflag, Bflag});

    // ARRAY_SHOW(A, L1, L2);
    // ARRAY_SHOW(B, L2, L1);

    int res   = 0;
    auto dhAm = get_ldims(Aflag, idx_vec).mapldims;
    auto dhBm = get_ldims(Bflag, idx_vec).mapldims;
    std::vector<int> idx(idx_vec.size(), 0);
    int imin = 0, imax = idx_vec.size() - 1;

    begin = std::chrono::steady_clock::now();
    {
        for (int loop = 0; loop < Nloop; ++loop) {
            res = 0;
            memset(idx.data(), 0, idx.size() * sizeof(int));
            int idxA = 0;
            int idxB = 0;
            while (true) {
                int i = imax;

                while (++idx[i] == idx_vec[i].d && i > imin) idx[i--] = 0;
                if (i == imin && idx[i] == idx_vec[i].d) break;

                idxA += dhAm[i];
                idxB += dhBm[i];

                res += A[idxA] * B[idxB];
            }
        }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "timing = " << static_cast<std::chrono::duration<double>>(end - begin).count() << "\n";
    std::cout << res << "!\n";

    begin = std::chrono::steady_clock::now();
    {
        EinsumHelper EH("cb,bc", {{L1, L2}, {L2, L1}});
        for (int loop = 0; loop < Nloop; ++loop) {
            res = 0;
            einsum(EH, {A, B}, &res);
        }
    }
    end = std::chrono::steady_clock::now();
    std::cout << "timing = " << static_cast<std::chrono::duration<double>>(end - begin).count() << "\n";
    std::cout << res << "&\n";


    auto EH = EinsumHelper("abcc[x],cbadf->fda", {{3, 4, 5, 5, 7}, {5, 4, 3, 8, 6}});

    return 0;
}
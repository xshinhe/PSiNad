#ifndef SPARSE_UTILS_H
#define SPARSE_UTILS_H
#include "definitions.h"

template <class T>
class SparseMat {
   public:
    enum { RowMajor, ColMajor };
    SparseMat(T* A, const int& N, const int& M, const int& major = RowMajor) {
        _major  = major;
        LDA     = (_major == RowMajor) ? N : M;
        _LLDA   = (_major == RowMajor) ? M : N;
        csr_LDA = new int[LDA];
        _nallow = 0;
        size    = 0;
        for (int i = 0; i < N; ++i) {  // Rowmajor !!!
            for (int j = 0; j < M; ++j) {
                if (NORM_OF(A[i * M + j]) > phys::math::eps8) { Push(i, j, A[i * M + j]); }
            }
        }
        Fix();
    }

    int Push(const int& n, const int& m, const T& val) {
        if (_major == RowMajor && n >= _nallow) {
            while (_nallow < n) {
                csr_LDA[_nallow] = size;
                _nallow++;
            }
            vec_val.push_back(val);
            vec_conn.push_back(m);
            size++;
            // std::cout << n << " " << m << ":" << val << std::endl;
        } else if (_major == ColMajor && m >= _nallow) {
            exit(-1);
        }
        return 0;
    }

    int Fix() {
        if (fixed) return -1;
        fixed    = true;
        csr_conn = new int[size];
        _val     = new T[size];
        for (int i = 0; i < size; ++i) {
            csr_conn[i] = vec_conn[i];
            _val[i]     = vec_val[i];
        }
        for (int n = _nallow; n < LDA; ++n) { csr_LDA[n] = size; }
        std::vector<int>().swap(vec_conn);
        std::vector<T>().swap(vec_val);
        return 0;
    }

    int MatMul(T* A, T* B) {
        plFunction();
        for (int i = 0; i < LDA; ++i) {
            A[i]   = (T) 0;
            int ib = (i == 0) ? 0 : csr_LDA[i - 1];
            int ie = csr_LDA[i];
            for (int ic = ib; ic < ie; ++ic) { A[i] += _val[ic] * B[csr_conn[ic]]; }
        }
        return 0;
    }

    int Show() {
        for (int i = 0; i < LDA; ++i) {
            int ib = (i == 0) ? 0 : csr_LDA[i - 1];
            int ie = csr_LDA[i];
            std::cout << std::setprecision(3) << std::setiosflags(std::ios::scientific) << setiosflags(std::ios::left)
                      << std::setw(7);
            int ic = ib;
            int j  = 0;
            while (ic < ie) {
                if (j < csr_conn[ic]) {
                    std::cout << "(X,X)\t";
                    j++;
                } else {
                    std::cout << _val[ic] << "\t";
                    ic++;
                    j++;
                }
            }
            while (j < _LLDA) {
                std::cout << "(X,X)\t";
                j++;
            }
            std::cout << std::endl;
        }
        return 0;
    }

    virtual ~SparseMat() {
        delete[] csr_LDA;
        if (fixed) { delete[] csr_conn, _val; }
    }

   private:
    bool fixed = false;
    int LDA, _LLDA, _major;
    int size;
    int* csr_LDA;
    int* csr_conn;
    int _nallow;
    T* _val;
    std::vector<int> vec_conn{};
    std::vector<T> vec_val{};
};


#endif  // SPARSE_UTILS_H
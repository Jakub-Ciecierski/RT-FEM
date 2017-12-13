#ifndef PROJECT_GPUMVSPARSEMULTIPLICATION_H
#define PROJECT_GPUMVSPARSEMULTIPLICATION_H

struct cusparseContext;
typedef struct cusparseContext *cusparseHandle_t;

struct cusparseMatDescr;
typedef struct cusparseMatDescr *cusparseMatDescr_t;

namespace rtfem {

template<class T>
class SparseMatrixCSR;

/**
 * Symmetry issues:
 * https://devtalk.nvidia.com/default/topic/506693/cusparse-matrix-vector-multiplication-algorithm/
 *
 * Given the fact that the transpose operation y=L^T*x
 * is 10x slower than non-transpose version y=L*x, the symmetric
 * property does not show up any performance gain
 *
 * @tparam T
 */
template<class T>
class GPUMVSparseMultiplication {
public:

    GPUMVSparseMultiplication();
    ~GPUMVSparseMultiplication();

    void PreSolve(const SparseMatrixCSR<T>& A);

    void Solve(T* x, T alpha,
               T* y, T beta);
private:
    void Terminate();

    int N;
    int nnz;
    int *d_col;
    int *d_row;
    T *d_val;
    T *d_x;
    T *d_y;

    cusparseHandle_t cusparseHandle;
    cusparseMatDescr_t description;
};
}


#endif //PROJECT_GPUMVSPARSEMULTIPLICATION_H

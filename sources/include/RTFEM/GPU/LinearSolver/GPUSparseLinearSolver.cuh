#ifndef PROJECT_GPUSPARSELINEARSOLVER_H
#define PROJECT_GPUSPARSELINEARSOLVER_H

struct cusparseContext;
typedef struct cusparseContext *cusparseHandle_t;

struct cublasContext;
typedef struct cublasContext *cublasHandle_t;

struct cusparseMatDescr;
typedef struct cusparseMatDescr *cusparseMatDescr_t;

namespace rtfem {

template<class T>
class SparseMatrixCSR;

template<class T>
class GPUSparseLinearSolver {
public:

    GPUSparseLinearSolver();

    virtual ~GPUSparseLinearSolver() = default;

    virtual void PreSolve(const SparseMatrixCSR<T>& A) = 0;

    virtual void Solve(const T* b, T* x) = 0;
protected:
    virtual void Terminate() = 0;

    int *d_col;
    int *d_row;
    T *d_val;
    int N;
    int nnz;

    cusparseHandle_t cusparseHandle;
    cublasHandle_t cublasHandle;
    cusparseMatDescr_t description;

    bool pre_solved_;
};
}


#endif //PROJECT_GPUSPARSELINEARSOLVER_H

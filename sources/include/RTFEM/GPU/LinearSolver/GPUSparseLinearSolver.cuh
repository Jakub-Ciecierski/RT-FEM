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
    virtual ~GPUSparseLinearSolver();

    virtual void PreSolve(const SparseMatrixCSR<T>& A);

    virtual void Solve(const T* b, T* x);
protected:
    virtual void Terminate();

    int *d_col;
    int *d_row;
    T *d_val;
    int N;
    int nnz;

    T *d_x;
    T *d_r;
    T *d_p;
    T *d_Ax;

    cusparseHandle_t cusparseHandle;
    cublasHandle_t cublasHandle;
    cusparseMatDescr_t description;

    bool pre_solved_;
};

}

#endif //PROJECT_GPUSPARSELINEARSOLVER_H

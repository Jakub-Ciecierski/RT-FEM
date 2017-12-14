#include "RTFEM/GPU/LinearSolver/GPUSparseCGLinearSolver.cuh"

#include <RTFEM/DataStructure/SparseMatrixCSR.h>

#include <stdlib.h>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>
#include <stdexcept>

namespace rtfem {

template<class T>
GPUSparseCGLinearSolver<T>::GPUSparseCGLinearSolver() :
        d_x(nullptr),
        d_r(nullptr),
        d_p(nullptr),
        d_Ax(nullptr){}

template<class T>
GPUSparseCGLinearSolver<T>::~GPUSparseCGLinearSolver(){
    Terminate();
}

template<class T>
void GPUSparseCGLinearSolver<T>::PreSolve(const SparseMatrixCSR<T>& A){
    this->pre_solved_ = true;
    this->N = A.n();
    this->nnz = A.values().size();

    cublasStatus_t cublasStatus;
    cudaError_t cuda_error;
    cublasStatus = cublasCreate(&this->cublasHandle);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&this->cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = cusparseCreateMatDescr(&this->description);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(this->description, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(this->description, CUSPARSE_INDEX_BASE_ZERO);

    cuda_error = cudaMalloc((void **)&this->d_col, this->nnz*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&this->d_row, (this->N+1)*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&this->d_val, this->nnz*sizeof(T));
    assert(cuda_error == cudaSuccess);

    cuda_error = cudaMalloc((void **)&d_x, this->N*sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_r, this->N*sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_p, this->N*sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_Ax, this->N*sizeof(T));
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(this->d_col, A.columns_indices().data(),
               this->nnz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_row, A.row_extents().data(),
               (this->N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_val, A.values().data(),
               this->nnz*sizeof(T), cudaMemcpyHostToDevice);
}

template<>
void GPUSparseCGLinearSolver<float>::PreSolve(const SparseMatrixCSR<float>& A){
    throw std::invalid_argument(
            "GPUSparseLinearSolver<float>::PreSolve not implemented");
}

template<class T>
void GPUSparseCGLinearSolver<T>::Solve(
        const T* B, T* x){
    const T tol = 1e-5f;
    const int max_iter = 1000;
    T a, b, na, r0, r1;
    T dot;
    int k;
    T alpha, beta, alpham1;

    cudaMemcpy(d_x, x, this->N*sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, B, this->N*sizeof(T), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.0;

    cusparseDcsrmv(this->cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                   this->N, this->N, this->nnz,
                   &alpha, this->description,
                   this->d_val, this->d_row, this->d_col, d_x, &beta, d_Ax);

    cublasDaxpy(this->cublasHandle, this->N,
                &alpham1, d_Ax, 1, d_r, 1);
    cublasDdot(this->cublasHandle,
               this->N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasDscal(this->cublasHandle, this->N, &b, d_p, 1);
            cublasDaxpy(this->cublasHandle, this->N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasDcopy(this->cublasHandle, this->N, d_r, 1, d_p, 1);
        }

        cusparseDcsrmv(this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                       this->N, this->N, this->nnz,
                       &alpha, this->description, this->d_val,
                       this->d_row, this->d_col, d_p,
                       &beta, d_Ax);
        cublasDdot(this->cublasHandle, this->N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasDaxpy(this->cublasHandle, this->N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasDaxpy(this->cublasHandle, this->N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasDdot(this->cublasHandle, this->N, d_r, 1, d_r, 1, &r1);
        cudaThreadSynchronize();
        k++;
    }

    cudaMemcpy(x, d_x, this->N*sizeof(T), cudaMemcpyDeviceToHost);
}

template<>
void GPUSparseCGLinearSolver<float>::Solve(
        const float* b, float* x){
    throw std::invalid_argument(
            "GPUSparseLinearSolver<float>::Solve not implemented");
}

template<class T>
void GPUSparseCGLinearSolver<T>::Terminate(){
    if(this->pre_solved_) {
        if(this->description)
            cusparseDestroyMatDescr(this->description);

        if(this->cusparseHandle)
            cusparseDestroy(this->cusparseHandle);
        if(this->cublasHandle)
            cublasDestroy(this->cublasHandle);
        if(this->d_col)
            cudaFree(this->d_col);
        if(this->d_row)
            cudaFree(this->d_row);
        if(this->d_val)
            cudaFree(this->d_val);
        if(d_x)
            cudaFree(d_x);
        if(d_r)
            cudaFree(d_r);
        if(d_p)
            cudaFree(d_p);
        if(d_Ax)
            cudaFree(d_Ax);
    }
}

template
class GPUSparseCGLinearSolver<float>;
template
class GPUSparseCGLinearSolver<double>;

}
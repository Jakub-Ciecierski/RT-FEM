#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

#include <RTFEM/DataStructure/SparseMatrixCSR.h>

#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>
#include <stdexcept>

namespace rtfem {

template<class T>
GPUSparseLinearSolver<T>::GPUSparseLinearSolver() :
        d_col(nullptr),
        d_row(nullptr),
        d_val(nullptr),
        N(0),
        nnz(0),
        d_x(nullptr),
        d_r(nullptr),
        d_p(nullptr),
        d_Ax(nullptr),
        cusparseHandle(nullptr),
        cublasHandle(nullptr),
        description(nullptr),
        pre_solved_(false){}

template<class T>
GPUSparseLinearSolver<T>::~GPUSparseLinearSolver(){
    Terminate();
}

template<class T>
void GPUSparseLinearSolver<T>::PreSolve(const SparseMatrixCSR<T>& A){
    pre_solved_ = true;
    N = A.n();
    nnz = A.values().size();

    cublasStatus_t cublasStatus;
    cudaError_t cuda_error;
    cublasStatus = cublasCreate(&cublasHandle);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = cusparseCreateMatDescr(&description);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(description, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(description, CUSPARSE_INDEX_BASE_ZERO);

    cuda_error = cudaMalloc((void **)&d_col, nnz*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_row, (N+1)*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_val, nnz*sizeof(double));
    assert(cuda_error == cudaSuccess);

    cuda_error = cudaMalloc((void **)&d_x, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_r, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_p, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_Ax, N*sizeof(double));
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(d_col, A.columns_indices().data(),
               nnz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, A.row_extents().data(),
               (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.values().data(),
               nnz*sizeof(double), cudaMemcpyHostToDevice);
}

template<>
void GPUSparseLinearSolver<float>::PreSolve(const SparseMatrixCSR<float>& A){
    throw std::invalid_argument(
            "GPUSparseLinearSolver<float>::PreSolve not implemented");
}

template<class T>
void GPUSparseLinearSolver<T>::Solve(
        const T* B, T* x){
    const T tol = 1e-5f;
    const int max_iter = 1000;
    T a, b, na, r0, r1;
    T dot;
    int k;
    T alpha, beta, alpham1;

    cudaMemcpy(d_x, x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, B, N*sizeof(double), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.0;

    cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                   N, N, nnz,
                   &alpha, description, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasDscal(cublasHandle, N, &b, d_p, 1);
            cublasDaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N,
                       N, nnz, &alpha, description, d_val, d_row, d_col, d_p,
                       &beta, d_Ax);
        cublasDdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasDaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasDaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaThreadSynchronize();
        k++;
    }

    cudaMemcpy(x, d_x, N*sizeof(double), cudaMemcpyDeviceToHost);
}

template<>
void GPUSparseLinearSolver<float>::Solve(
        const float* b, float* x){
    throw std::invalid_argument(
            "GPUSparseLinearSolver<float>::Solve not implemented");
}

template<class T>
void GPUSparseLinearSolver<T>::Terminate(){
    if(pre_solved_) {
        if(cusparseHandle)
            cusparseDestroy(cusparseHandle);
        if(cublasHandle)
            cublasDestroy(cublasHandle);
        if(d_col)
            cudaFree(d_col);
        if(d_row)
            cudaFree(d_row);
        if(d_val)
            cudaFree(d_val);
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
class GPUSparseLinearSolver<float>;
template
class GPUSparseLinearSolver<double>;

}
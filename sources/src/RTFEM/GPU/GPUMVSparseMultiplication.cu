#include "RTFEM/GPU/GPUMVSparseMultiplication.cuh"

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>

#include <RTFEM/DataStructure/SparseMatrixCSR.h>

namespace rtfem {

template<class T>
void GPUMVSparseMultiplication<T>::Solve(
        const SparseMatrixCSR<T>& A,
        T* x, T alpha,
        T* y, T beta){
    int *d_col, *d_row;
    double *d_val;
    double *d_x;
    double *d_y;

    int nz = A.values().size();
    int N = A.n();

    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cudaError_t cuda_error;

    cublasStatus = cublasCreate(&cublasHandle);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;

    cusparseStatus = cusparseCreate(&cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    cuda_error = cudaMalloc((void **)&d_col, nz*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_row, (N+1)*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_val, nz*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&y, N*sizeof(double));
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(d_col, A.columns_indices().data(),
               nz*sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, A.row_extents().data(),
               (N+1)*sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.values().data(),
               nz*sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, N*sizeof(double), cudaMemcpyHostToDevice);

    cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                   N, N, nz,
                   &alpha,
                   descr,
                   d_val, d_row, d_col,
                   d_x,
                   &beta, d_y);

    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_y);
}

}
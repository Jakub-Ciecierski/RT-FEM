#include "RTFEM/GPU/GPUMVSparseMultiplication.cuh"

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>

#include <RTFEM/DataStructure/SparseMatrixCSR.h>
#include <stdexcept>

namespace rtfem {

template<class T>
GPUMVSparseMultiplication<T>::GPUMVSparseMultiplication() : N(0),
                                                            nnz(0),
                                                            d_col(nullptr),
                                                            d_row(nullptr),
                                                            d_val(nullptr),
                                                            d_x(nullptr),
                                                            d_y(nullptr) {}

template<class T>
GPUMVSparseMultiplication<T>::~GPUMVSparseMultiplication(){
    Terminate();
}

template<class T>
void GPUMVSparseMultiplication<T>::PreSolve(const SparseMatrixCSR<T>& A){
    nnz = A.values().size();
    N = A.n();

    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = cusparseCreateMatDescr(&description);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    switch(A.type()){
        case MatrixType::General:
            cusparseSetMatType(description, CUSPARSE_MATRIX_TYPE_GENERAL);
            break;
        case MatrixType::Symmetric:
            cusparseSetMatType(description, CUSPARSE_MATRIX_TYPE_SYMMETRIC);
            break;
    }
    cusparseSetMatIndexBase(description, CUSPARSE_INDEX_BASE_ZERO);

    cudaError_t cuda_error;
    cuda_error = cudaMalloc((void **)&d_col, nnz*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_row, (N+1)*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_val, nnz*sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, N*sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_y, N*sizeof(T));
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(d_col, A.columns_indices().data(),
               nnz*sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, A.row_extents().data(),
               (N+1)*sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, A.values().data(),
               nnz*sizeof(T),
               cudaMemcpyHostToDevice);
}

template<>
void GPUMVSparseMultiplication<float>::PreSolve(
        const SparseMatrixCSR<float>& A){
    throw std::invalid_argument(
            "GPUMVSparseMultiplication<float>::Solve not implemented");
}

template<class T>
void GPUMVSparseMultiplication<T>::Solve(
        T* x, T alpha,
        T* y, T beta){
    cudaError_t cuda_error;

    cuda_error = cudaMemcpy(d_x, x, N*sizeof(T), cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemcpy(d_y, y, N*sizeof(T), cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);

    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseDcsrmv(cusparseHandle,
                                    CUSPARSE_OPERATION_NON_TRANSPOSE,
                                    N, N, nnz,
                                    &alpha,
                                    description,
                                    d_val, d_row, d_col,
                                    d_x,
                                    &beta,
                                    d_y);
    cuda_error = cudaDeviceSynchronize();
    assert(cuda_error == cudaSuccess);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cuda_error = cudaMemcpy(y, d_y, sizeof(T)*N,
                            cudaMemcpyDeviceToHost);
    assert(cuda_error == cudaSuccess);
}

template<>
void GPUMVSparseMultiplication<float>::Solve(
        float* x, float alpha,
        float* y, float beta){
    throw std::invalid_argument(
            "GPUMVSparseMultiplication<float>::Solve not implemented");
}

template<class T>
void GPUMVSparseMultiplication<T>::Terminate(){
    cusparseDestroy(cusparseHandle);

    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_y);
}

template
class GPUMVSparseMultiplication<double>;
template
class GPUMVSparseMultiplication<float>;

}
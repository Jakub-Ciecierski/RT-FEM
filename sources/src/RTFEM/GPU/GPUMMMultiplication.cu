#include "RTFEM/GPU/GPUMMMultiplication.cuh"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <assert.h>
#include <stdexcept>

namespace rtfem {

template <class T>
GPUMMMultiplication<T>::GPUMMMultiplication(){}

template <class T>
GPUMMMultiplication<T>::~GPUMMMultiplication(){}

template <class T>
void GPUMMMultiplication<T>::Solve(const T* A, const T* B, T* C,
                                   T alpha, T beta,
                                   int m, int k, int n,
                                   MatrixOperation A_operation,
                                   MatrixOperation B_operation){
    cudaError_t cuda_error;
    cublasStatus_t status;
    cublasHandle_t handle;

    auto GetOperation = [](const MatrixOperation& operation){
        switch(operation){
            case MatrixOperation::None:
                return CUBLAS_OP_N;
            case MatrixOperation::Transpose:
                return CUBLAS_OP_T;
            default:
                return CUBLAS_OP_N;
        }
    };

    T *d_A = nullptr;
    T *d_B = nullptr;
    T *d_C = nullptr;

    cuda_error = cudaMalloc((void **) &d_A, m * k * sizeof(*A));
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc((void **) &d_B, k * n * sizeof(*B));
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc((void **) &d_C, m * n * sizeof(*C));
    assert(cudaSuccess == cuda_error);

    status = cublasCreate(&handle);
    assert(CUBLAS_STATUS_SUCCESS == status);
    status = cublasSetMatrix(m, k, sizeof(*A), A, m, d_A, m);
    assert(CUBLAS_STATUS_SUCCESS == status);
    status = cublasSetMatrix(k, n, sizeof(*B), B, k, d_B, k);
    assert(CUBLAS_STATUS_SUCCESS == status);
    status = cublasSetMatrix(m, n, sizeof(*C), C, m, d_C, m);
    assert(CUBLAS_STATUS_SUCCESS == status);

    status = cublasDgemm(handle,
                         GetOperation(A_operation),
                         GetOperation(B_operation),
                         m, n, k,
                         &alpha,
                         d_A, m, d_B, k,
                         &beta, d_C, m);
    assert(CUBLAS_STATUS_SUCCESS == status);

    status = cublasGetMatrix(m, n, sizeof(*C), d_C, m, C, m);
    assert(CUBLAS_STATUS_SUCCESS == status);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    cublasDestroy(handle);
}

template<>
void GPUMMMultiplication<float>::Solve(const float *A, const float *B, float *C,
                                       float alpha, float beta,
                                       int m, int k, int n,
                                       MatrixOperation A_operation,
                                       MatrixOperation B_operation) {
    throw std::invalid_argument(
            "GPUMMMultiplication<float>::Solve not implemented");
}

template
class GPUMMMultiplication<double>;
template
class GPUMMMultiplication<float>;

}
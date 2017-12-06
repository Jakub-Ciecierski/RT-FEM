#include "RTFEM/GPU/GPUMatrixMultiplication.cuh"

#include <cublas_v2.h>
#include <cuda_runtime.h>

#include <cstdlib>
#include <assert.h>
#include <stdexcept>

namespace rtfem {

template<class T>
GPUMatrixMultiplication<T>::GPUMatrixMultiplication() : d_A_(nullptr), n_(0) {}

template<class T>
GPUMatrixMultiplication<T>::~GPUMatrixMultiplication(){
    Terminate();
}

template<class T>
void GPUMatrixMultiplication<T>::PreSolve(T* A, int n){
    n_ = n;
    cudaError_t cuda_error;
    cublasStatus_t status;

    cuda_error = cudaMalloc((void **) &d_A_, n_ * n_ * sizeof(*A));
    assert(cudaSuccess == cuda_error);

    status = cublasCreate(&cublas_handle_);
    assert(CUBLAS_STATUS_SUCCESS == status);

    status = cublasSetMatrix(n_, n_, sizeof(*A), A, n_, d_A_, n_);
    assert(CUBLAS_STATUS_SUCCESS == status);
}

template<>
void GPUMatrixMultiplication<float>::PreSolve(float* A, int n){
    throw std::invalid_argument(
            "GPUMatrixMultiplication<float>::PreSolve not implemented");
}

template<class T>
void GPUMatrixMultiplication<T>::Solve(T* x, T alpha,
                                       T* y, T beta){
    T *d_x = nullptr;
    T *d_y = nullptr;

    cudaError_t cuda_error;
    cublasStatus_t status;

    cuda_error = cudaMalloc((void **) &d_x, n_ * sizeof(*x));
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc((void **) &d_y, n_ * sizeof(*y));
    assert(cudaSuccess == cuda_error);

    status = cublasSetVector(n_, sizeof(*x), x, 1, d_x, 1);
    assert(CUBLAS_STATUS_SUCCESS == status);
    status = cublasSetVector(n_, sizeof(*y), y, 1, d_y, 1);
    assert(CUBLAS_STATUS_SUCCESS == status);

    status = cublasDgemv(cublas_handle_, CUBLAS_OP_N,
                         n_, n_,
                         &alpha, d_A_, n_, d_x, 1,
                         &beta, d_y, 1);
    assert(CUBLAS_STATUS_SUCCESS == status);

    status = cublasGetVector(n_, sizeof(*y), d_y, 1, y, 1);
    assert(CUBLAS_STATUS_SUCCESS == status);

    cudaFree(d_x);
    cudaFree(d_y);
}

template<>
void GPUMatrixMultiplication<float>::Solve(float* x, float alpha,
                                           float* y, float beta){
    throw std::invalid_argument(
            "GPUMatrixMultiplication<float>::Solve not implemented");
}

template<class T>
void GPUMatrixMultiplication<T>::Terminate(){
    cudaFree(d_A_);
    cublasDestroy(cublas_handle_);
}

template
class GPUMatrixMultiplication<double>;
template
class GPUMatrixMultiplication<float>;

}

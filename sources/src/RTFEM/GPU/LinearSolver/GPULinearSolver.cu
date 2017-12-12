#include "RTFEM/GPU/LinearSolver/GPULinearSolver.cuh"

#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>

#include <cusolverDn.h>

namespace rtfem {

template<class T>
GPULinearSolver<T>::GPULinearSolver() :
    d_A(nullptr),
    d_b(nullptr),
    d_pivot(nullptr),
    d_info(nullptr),
    d_work(nullptr),
    cusolverH(nullptr),
    stream(nullptr) {}

template <class T>
GPULinearSolver<T>::~GPULinearSolver(){
    Terminate();
}

template <class T>
CUDA_HOST_MEMBER
void GPULinearSolver<T>::PreSolve(T* A, int n){
    n_ = n;

    // Host
    int info = 0;
    int lwork = 0;

    cudaError_t cuda_error = cudaSuccess;

    /* step 1: create cusolver handle, bind a stream */
    auto status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cuda_error = cudaStreamCreateWithFlags(&stream,
                                           cudaStreamNonBlocking);
    assert(cudaSuccess == cuda_error);

    status = cusolverDnSetStream(cusolverH, stream);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    /* step 2: copy A to device */
    cuda_error = cudaMalloc ((void**)&d_A, sizeof(T) * n * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_b, sizeof(T) * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_pivot, sizeof(int) * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_info, sizeof(int));
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(d_A, A, sizeof(T)*n*n, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cuda_error);

    /* step 3: query working space of getrf */
    status = cusolverDnDgetrf_bufferSize(
        cusolverH,
        n,
        n,
        d_A,
        n,
        &lwork);
    assert(CUSOLVER_STATUS_SUCCESS == status);

    cuda_error = cudaMalloc((void**)&d_work, sizeof(T)*lwork);
    assert(cudaSuccess == cuda_error);

    /* step 4: LU factorization */
    status = cusolverDnDgetrf(
        cusolverH,
        n,
        n,
        d_A,
        n,
        d_work,
        d_pivot,
        d_info);
    cuda_error = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);

    if ( 0 > info ){
        printf("%d-th parameter is wrong \n", -info);
        exit(1);
    }

}

template <>
CUDA_HOST_MEMBER
void GPULinearSolver<float>::PreSolve(float* A, int n){
    throw std::invalid_argument(
        "GPULinearSolver<float>::PreSolve not implemented");
}

template <class T>
CUDA_HOST_MEMBER
void GPULinearSolver<T>::Solve(const T* b, int n, T* x){
    cudaError_t cuda_error = cudaSuccess;

    cuda_error = cudaMemcpy(d_b, b, sizeof(T)*n, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cuda_error);

    auto status = cusolverDnDgetrs(
        cusolverH,
        CUBLAS_OP_N,
        n,
        1, /* nrhs */
        d_A,
        n,
        d_pivot,
        d_b,
        n,
        d_info);

    cuda_error = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == status);
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(x, d_b, sizeof(T)*n,
                            cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);
}

template <>
CUDA_HOST_MEMBER
void GPULinearSolver<float>::Solve(const float* b, int n, float* x){
    throw std::invalid_argument(
            "GPULinearSolver<float>::Solve not implemented");
}

template <class T>
CUDA_HOST_MEMBER
void GPULinearSolver<T>::Terminate(){
    if (d_A    ) cudaFree(d_A);
    if (d_b    ) cudaFree(d_b);
    if (d_pivot) cudaFree(d_pivot);
    if (d_info ) cudaFree(d_info);
    if (d_work ) cudaFree(d_work);

    if (cusolverH   ) cusolverDnDestroy(cusolverH);
    if (stream      ) cudaStreamDestroy(stream);
}

template
class GPULinearSolver<double>;
template
class GPULinearSolver<float>;

}
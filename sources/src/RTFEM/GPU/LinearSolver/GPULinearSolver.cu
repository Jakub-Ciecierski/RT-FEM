#include "RTFEM/GPU/LinearSolver/GPULinearSolver.cuh"

#include <cusolverDn.h>

#include <assert.h>
#include <cstdlib>
#include <cstdio>

namespace rtfem {

CUDA_HOST_MEMBER
void GPULinearSolver::Solve(double* A, double* b, int n, double* x){
    cusolverDnHandle_t cusolverH = nullptr;
    cudaStream_t stream = nullptr;

    // Host
    double* LU = (double*)malloc(sizeof(double) * n * n);
    int* Ipivot = (int*)malloc(sizeof(int) * n);
    int info = 0;

    // Device
    double* d_A = nullptr;
    double* d_b = nullptr;
    int* d_pivot = nullptr;
    int* d_info = nullptr;
    int lwork = 0;
    double* d_work = nullptr;
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
    cuda_error = cudaMalloc ((void**)&d_A, sizeof(double) * n * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_b, sizeof(double) * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_pivot, sizeof(int) * n);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMalloc ((void**)&d_info, sizeof(int));
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(d_A, A, sizeof(double)*n*n, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cuda_error);
    cuda_error = cudaMemcpy(d_b, b, sizeof(double)*n, cudaMemcpyHostToDevice);
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

    cuda_error = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
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

    cuda_error = cudaMemcpy(Ipivot , d_pivot, sizeof(int)*n,
                            cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(LU, d_A, sizeof(double)*n*n,
                           cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);

    cuda_error = cudaMemcpy(&info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);

    if ( 0 > info ){
        printf("%d-th parameter is wrong \n", -info);
        exit(1);
    }


    /*
     * step 5: solve A*X = B
     *       | 1 |       | -0.3333 |
     *   B = | 2 |,  X = |  0.6667 |
     *       | 3 |       |  0      |
     *
     */
    status = cusolverDnDgetrs(
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

    cuda_error = cudaMemcpy(x, d_b, sizeof(double)*n,
                            cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cuda_error);

    /* free resources */
    if (d_A    ) cudaFree(d_A);
    if (d_b    ) cudaFree(d_b);
    if (d_pivot) cudaFree(d_pivot);
    if (d_info ) cudaFree(d_info);
    if (d_work ) cudaFree(d_work);

    if(LU)
        free(LU);
    if(Ipivot)
        free(Ipivot);

    if (cusolverH   ) cusolverDnDestroy(cusolverH);
    if (stream      ) cudaStreamDestroy(stream);

    //cudaDeviceReset();
}


}
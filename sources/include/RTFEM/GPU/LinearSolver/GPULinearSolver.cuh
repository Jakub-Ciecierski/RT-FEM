#ifndef PROJECT_LINEARSOLVER_H
#define PROJECT_LINEARSOLVER_H

#include <RTFEM/GPU/CudaFunctionMacros.cuh>

struct cusolverDnContext;
typedef struct cusolverDnContext *cusolverDnHandle_t;

struct CUstream_st;
typedef __device_builtin__ struct CUstream_st *cudaStream_t;

namespace rtfem {

/**
 * Uses cusolver library
 * @tparam T
 */
template <class T>
class GPULinearSolver {
public:
    GPULinearSolver();
    ~GPULinearSolver();

    CUDA_HOST_MEMBER void PreSolve(T* A, int n);
    CUDA_HOST_MEMBER void Solve(const T* b, int n, T* x);

private:
    CUDA_HOST_MEMBER void Terminate();

    T* d_A;
    T* d_b;

    int* d_pivot;
    int* d_info;

    T* d_work = nullptr;

    int n_;

    cusolverDnHandle_t cusolverH;
    cudaStream_t stream;
};

}

#endif //PROJECT_LINEARSOLVER_H

#ifndef PROJECT_GPUMATRIXMULTIPLICATION_H
#define PROJECT_GPUMATRIXMULTIPLICATION_H

struct cublasContext;
typedef struct cublasContext *cublasHandle_t;

namespace rtfem {

/**
 * Uses cublas library
 */
template<class T>
class GPUMVMultiplication {
public:

    GPUMVMultiplication();
    ~GPUMVMultiplication();

    void PreSolve(T* A, int n);

    void Solve(T *x, T scalar,
               T *y, T beta);
private:
    void Terminate();

    T* d_A_;
    int n_;

    cublasHandle_t cublas_handle_;
};
}


#endif //PROJECT_GPUMATRIXMULTIPLICATION_H

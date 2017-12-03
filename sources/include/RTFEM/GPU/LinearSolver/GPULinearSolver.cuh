#ifndef PROJECT_LINEARSOLVER_H
#define PROJECT_LINEARSOLVER_H

#include <RTFEM/GPU/CudaFunctionMacros.cuh>

namespace rtfem {

template <class T>
class GPULinearSolver {
public:
    GPULinearSolver() = default;
    ~GPULinearSolver() = default;

    CUDA_HOST_MEMBER void Solve(T* A, T* b, int n, T* x);
private:
};

}

#endif //PROJECT_LINEARSOLVER_H

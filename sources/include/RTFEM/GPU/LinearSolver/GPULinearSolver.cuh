#ifndef PROJECT_LINEARSOLVER_H
#define PROJECT_LINEARSOLVER_H

#include <RTFEM/GPU/CudaFunctionMacros.cuh>

namespace rtfem {

class GPULinearSolver {
public:
    GPULinearSolver() = default;
    ~GPULinearSolver() = default;

    CUDA_HOST_MEMBER void Solve(double* A, double* b, int n, double* x);
private:
};

}

#endif //PROJECT_LINEARSOLVER_H

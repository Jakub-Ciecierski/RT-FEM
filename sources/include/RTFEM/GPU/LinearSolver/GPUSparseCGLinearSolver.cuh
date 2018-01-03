#ifndef PROJECT_GPUSPARSECGLINEARSOLVER_H
#define PROJECT_GPUSPARSECGLINEARSOLVER_H

#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

namespace rtfem {

/**
 * Cuda sample
 * @tparam T
 */
template<class T>
class GPUSparseCGLinearSolver : public GPUSparseLinearSolver<T>{
public:
    GPUSparseCGLinearSolver();
    virtual ~GPUSparseCGLinearSolver();

    virtual void PreSolve(const SparseMatrixCSR<T>& A) override;

    virtual void Solve(const T* b, T* x) override;

    void SetTolerance(T tolerance);
protected:
    virtual void Terminate() override;

    T *d_x;
    T *d_r;
    T *d_p;
    T *d_Ax;

    T tolerance_;
};

}

#endif //PROJECT_GPUSPARSECGLINEARSOLVER_H

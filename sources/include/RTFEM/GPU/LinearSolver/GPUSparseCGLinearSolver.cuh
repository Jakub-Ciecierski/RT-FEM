#ifndef PROJECT_GPUSPARSECGLINEARSOLVER_H
#define PROJECT_GPUSPARSECGLINEARSOLVER_H

#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

namespace rtfem {

template<class T>
class GPUSparseCGLinearSolver : public GPUSparseLinearSolver<T>{
public:
    GPUSparseCGLinearSolver();
    virtual ~GPUSparseCGLinearSolver();

    virtual void PreSolve(const SparseMatrixCSR<T>& A) override;

    virtual void Solve(const T* b, T* x) override;
protected:
    virtual void Terminate() override;

    T *d_x;
    T *d_r;
    T *d_p;
    T *d_Ax;

};

}

#endif //PROJECT_GPUSPARSECGLINEARSOLVER_H

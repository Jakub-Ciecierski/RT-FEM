#ifndef PROJECT_GPUSPARSEPRECONDLINEARSOLVER_H
#define PROJECT_GPUSPARSEPRECONDLINEARSOLVER_H

#include "RTFEM/GPU/LinearSolver/GPUSparseCGLinearSolver.cuh"

struct cusparseSolveAnalysisInfo;
typedef struct cusparseSolveAnalysisInfo *cusparseSolveAnalysisInfo_t;

struct cusparseSolveAnalysisInfo;
typedef struct cusparseSolveAnalysisInfo *cusparseSolveAnalysisInfo_t;

namespace rtfem {

/**
 * Cuda sample
 * @tparam T
 */
template<class T>
class GPUSparseCGPrecondLinearSolver : public GPUSparseCGLinearSolver<T>{
public:

    GPUSparseCGPrecondLinearSolver();
    ~GPUSparseCGPrecondLinearSolver();

    virtual void PreSolve(const SparseMatrixCSR<T>& A) override;

    virtual void Solve(const T* b, T* x) override;
protected:
    virtual void Terminate() override;

private:
    T* d_y;
    T* d_omega;
    T* d_valsILU0;
    T* d_zm1;
    T* d_zm2;
    T* d_rm2;

    cusparseSolveAnalysisInfo_t infoA;
    cusparseSolveAnalysisInfo_t info_u;
    cusparseMatDescr_t descrL;
    cusparseMatDescr_t descrU;
};
}


#endif //PROJECT_GPUSPARSEPRECONDLINEARSOLVER_H

#ifndef PROJECT_GPUSPARSELUSOLVER_H
#define PROJECT_GPUSPARSELUSOLVER_H

#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

struct csrsv2Info;
typedef struct csrsv2Info *csrsv2Info_t;

namespace rtfem {

/**
 * https://github.com/OrangeOwlSolutions/Linear-Algebra/blob/master/LU/LU_linear_system_sparse.cu
 * @tparam T
 */
template<class T>
class GPUSparseLULinearSolver : public GPUSparseLinearSolver<T>{
public:

    GPUSparseLULinearSolver();
    ~GPUSparseLULinearSolver();

    virtual void PreSolve(const SparseMatrixCSR<T>& A) override;

    virtual void Solve(const T* b, T* x) override;
protected:
    virtual void Terminate() override;

private:
    T *d_x;
    T *d_z;
    T *d_y;

    void *pBuffer = 0;

    cusparseMatDescr_t descr_L;
    cusparseMatDescr_t descr_U;

    csrsv2Info_t info_L;
    csrsv2Info_t info_U;
};
}


#endif //PROJECT_GPUSPARSELUSOLVER_H

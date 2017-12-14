#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

namespace rtfem {

template<class T>
GPUSparseLinearSolver<T>::GPUSparseLinearSolver() :
        d_col(nullptr),
        d_row(nullptr),
        d_val(nullptr),
        N(0),
        nnz(0),
        cusparseHandle(nullptr),
        cublasHandle(nullptr),
        description(nullptr),
        pre_solved_(false) {}

template
class GPUSparseLinearSolver<float>;
template
class GPUSparseLinearSolver<double>;

}
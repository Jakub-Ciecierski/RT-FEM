#include "RTFEM/GPU/LinearSolver/GPUSparseCGPrecondLinearSolver.cuh"

#include <RTFEM/DataStructure/SparseMatrixCSR.h>

#include <stdlib.h>
#include <stdio.h>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>
#include <stdexcept>
#include <iostream>

namespace rtfem {

template<class T>
GPUSparseCGPrecondLinearSolver<T>::GPUSparseCGPrecondLinearSolver() :
        d_y(nullptr),
        d_omega(nullptr),
        d_valsILU0(nullptr),
        d_zm1(nullptr),
        d_zm2(nullptr),
        d_rm2(nullptr),
        infoA(nullptr),
        info_u(nullptr),
        descrL(nullptr),
        descrU(nullptr) {}

template<class T>
GPUSparseCGPrecondLinearSolver<T>::~GPUSparseCGPrecondLinearSolver(){
    Terminate();
}

template<class T>
void GPUSparseCGPrecondLinearSolver<T>::PreSolve(const SparseMatrixCSR<T>& A){
    this->pre_solved_ = true;

    this->N = A.n();
    this->nnz = A.values().size();
    
    cublasStatus_t cublasStatus;
    cublasStatus = cublasCreate(&this->cublasHandle);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    /* Create CUSPARSE context */
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&this->cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /* Description of the A matrix*/
    cusparseStatus = cusparseCreateMatDescr(&this->description);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /* Define the properties of the matrix */
    cusparseSetMatType(this->description, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(this->description, CUSPARSE_INDEX_BASE_ZERO);

    /* Allocate required memory */
    cudaMalloc((void **)&this->d_col, this->nnz*sizeof(int));
    cudaMalloc((void **)&this->d_row, (this->N+1)*sizeof(int));
    cudaMalloc((void **)&this->d_val, this->nnz*sizeof(T));
    cudaMalloc((void **)&this->d_x, this->N*sizeof(T));
    cudaMalloc((void **)&this->d_y, this->N*sizeof(T));
    cudaMalloc((void **)&this->d_r, this->N*sizeof(T));
    cudaMalloc((void **)&this->d_p, this->N*sizeof(T));
    cudaMalloc((void **)&this->d_omega, this->N*sizeof(T));
    cudaMalloc((void **)&this->d_valsILU0, this->nnz*sizeof(T));
    cudaMalloc((void **)&this->d_zm1, (this->N)*sizeof(T));
    cudaMalloc((void **)&this->d_zm2, (this->N)*sizeof(T));
    cudaMalloc((void **)&this->d_rm2, (this->N)*sizeof(T));
    
    cudaMemcpy(this->d_col, A.columns_indices().data(), this->nnz * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_row, A.row_extents().data(), (this->N + 1) * sizeof(int),
               cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_val, A.values().data(), this->nnz * sizeof(T),
               cudaMemcpyHostToDevice);
    
    /* create the analysis info object for t`he A matrix */
    cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /* Perform the analysis for the Non-Transpose case */
    cusparseStatus = cusparseDcsrsv_analysis(this->cusparseHandle,
                                             CUSPARSE_OPERATION_NON_TRANSPOSE,
                                             this->N,
                                             this->nnz,
                                             this->description,
                                             this->d_val,
                                             this->d_row,
                                             this->d_col,
                                             infoA);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /* Copy A data to ILU0 vals as input*/
    cudaMemcpy(d_valsILU0,
               this->d_val,
               this->nnz*sizeof(T), cudaMemcpyDeviceToDevice);

    /* generate the Incomplete LU factor H for the matrix A using cudsparseScsrilu0 */
    cusparseStatus = cusparseDcsrilu0(this->cusparseHandle,
                                      CUSPARSE_OPERATION_NON_TRANSPOSE,
                                      this->N,
                                      this->description,
                                      d_valsILU0,
                                      this->d_row,
                                      this->d_col,
                                      infoA);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /* Create info objects for the ILU0 preconditioner */
    cusparseCreateSolveAnalysisInfo(&info_u);

    cusparseStatus = cusparseCreateMatDescr(&descrL);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);
    cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT);

    cusparseStatus = cusparseCreateMatDescr(&descrU);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(descrU,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cusparseStatus = cusparseDcsrsv_analysis(this->cusparseHandle,
                                             CUSPARSE_OPERATION_NON_TRANSPOSE,
                                             this->N,
                                             this->nnz,
                                             descrU,
                                             this->d_val,
                                             this->d_row,
                                             this->d_col,
                                             info_u);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);
}

template<>
void GPUSparseCGPrecondLinearSolver<float>::PreSolve(
        const SparseMatrixCSR<float>& A){
    throw std::invalid_argument(
            "GPUSparsePreCondLinearSolver<float>::PreSolve not implemented");
}

template<class T>
void GPUSparseCGPrecondLinearSolver<T>::Solve(const T* b, T* x){
    const int max_iter = 1000;
    const T tol = 1e-5f;

    int k = 0;
    T r1, alpha, beta;
    T numerator, denominator, nalpha;
    const T floatone = 1.0;
    const T floatzero = 0.0;

    cudaMemcpy(this->d_x, x, this->N*sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_r, b, this->N*sizeof(T), cudaMemcpyHostToDevice);

    int nzILU0 = 2*this->N-1;

    for (int i = 0; i < this->N; i++)
    {
        x[i] = 0.0;
    }

    cudaMemcpy(this->d_r, b, this->N*sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_x, x, this->N*sizeof(T), cudaMemcpyHostToDevice);

    k = 0;
    cublasDdot(this->cublasHandle,
               this->N,
               this->d_r, 1,
               this->d_r, 1, &r1);
    cusparseStatus_t cusparseStatus;
    while (r1 > tol*tol && k <= max_iter)
    {
        // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
        cusparseStatus = cusparseDcsrsv_solve(this->cusparseHandle,
                                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                                              this->N, &floatone, descrL,
                                              d_valsILU0,
                                              this->d_row,
                                              this->d_col,
                                              this->infoA,
                                              this->d_r, d_y);
        assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

        // Back Substitution
        cusparseStatus = cusparseDcsrsv_solve(this->cusparseHandle,
                                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                                              this->N, &floatone, descrU,
                                              d_valsILU0,
                                              this->d_row,
                                              this->d_col,
                                              this->info_u,
                                              this->d_y, d_zm1);
        assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

        k++;

        if (k == 1)
        {
            cublasDcopy(this->cublasHandle, this->N, d_zm1, 1,
                        this->d_p, 1);
        }
        else
        {
            cublasDdot(this->cublasHandle,
                       this->N,
                       this->d_r, 1, d_zm1, 1, &numerator);
            cublasDdot(this->cublasHandle,
                       this->N, d_rm2, 1, d_zm2, 1, &denominator);
            beta = numerator/denominator;
            cublasDscal(this->cublasHandle, this->N, &beta, this->d_p, 1);
            cublasDaxpy(this->cublasHandle, this->N,
                        &floatone, d_zm1, 1, this->d_p, 1) ;
        }

        cusparseDcsrmv(this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                       this->N, this->N,
                       nzILU0, &floatone, descrU,
                       this->d_val, this->d_row, this->d_col, this->d_p,
                       &floatzero, d_omega);
        cublasDdot(this->cublasHandle,
                   this->N, this->d_r, 1, d_zm1, 1, &numerator);
        cublasDdot(this->cublasHandle,
                   this->N, this->d_p, 1, d_omega, 1, &denominator);
        alpha = numerator / denominator;
        cublasDaxpy(this->cublasHandle,
                    this->N, &alpha,
                    this->d_p, 1, this->d_x, 1);
        cublasDcopy(this->cublasHandle,
                    this->N, this->d_r, 1, d_rm2, 1);
        cublasDcopy(this->cublasHandle,
                    this->N, d_zm1, 1, d_zm2, 1);
        nalpha = -alpha;
        cublasDaxpy(this->cublasHandle, this->N, &nalpha,
                    d_omega, 1, this->d_r, 1);
        cublasDdot(this->cublasHandle, this->N, this->d_r, 1,
                   this->d_r, 1, &r1);
    }

    cudaMemcpy(x, this->d_x, this->N*sizeof(T), cudaMemcpyDeviceToHost);
}

template<>
void GPUSparseCGPrecondLinearSolver<float>::Solve(
        const float* B, float* x){
    throw std::invalid_argument(
            "GPUSparsePreCondLinearSolver<float>::Solve not implemented");
}

template<class T>
void GPUSparseCGPrecondLinearSolver<T>::Terminate(){
    if(this->pre_solved_) {
        /* Destroy parameters */
        cusparseDestroySolveAnalysisInfo(infoA);
        cusparseDestroySolveAnalysisInfo(info_u);

        cusparseDestroyMatDescr(descrL);
        cusparseDestroyMatDescr(descrU);

        if(this->d_y)
            cudaFree(this->d_y);
        if(this->d_omega)
            cudaFree(this->d_omega);
        if(this->d_valsILU0)
            cudaFree(this->d_valsILU0);
        if(this->d_zm1)
            cudaFree(this->d_zm1);
        if(this->d_zm2)
            cudaFree(this->d_zm2);
        if(this->d_rm2)
            cudaFree(this->d_rm2);
    }
}

template
class GPUSparseCGPrecondLinearSolver<double>;
template
class GPUSparseCGPrecondLinearSolver<float>;

}
#include "RTFEM/GPU/LinearSolver/GPUSparseLULinearSolver.cuh"

#include <RTFEM/DataStructure/SparseMatrixCSR.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <math_functions.h>

namespace rtfem {

template<class T>
GPUSparseLULinearSolver<T>::GPUSparseLULinearSolver() :
        d_x(nullptr),
        d_z(nullptr),
        d_y(nullptr),
        pBuffer(nullptr),
        descr_L(nullptr),
        descr_U(nullptr),
        info_L(nullptr),
        info_U(nullptr){}

template<class T>
GPUSparseLULinearSolver<T>::~GPUSparseLULinearSolver(){}

template<class T>
void GPUSparseLULinearSolver<T>::PreSolve(const SparseMatrixCSR<T>& A){
    this->pre_solved_ = true;

    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&this->cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /*******************************/
    /* FROM DENSE TO SPARSE MATRIX */
    /*******************************/
    // --- Descriptor for sparse matrix A
    cusparseStatus = cusparseCreateMatDescr(&this->description);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);
    cusparseStatus = cusparseSetMatType(this->description,
                                        CUSPARSE_MATRIX_TYPE_GENERAL);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);
    cusparseStatus = cusparseSetMatIndexBase(this->description,
                                             CUSPARSE_INDEX_BASE_ZERO);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    this->N = A.n();
    this->nnz = A.values().size();

    cudaError_t cuda_error;
    // --- Device side sparse matrix
    cuda_error = cudaMalloc(&this->d_val, this->nnz * sizeof(T));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc(&this->d_row, (this->N + 1) * sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc(&this->d_col, this->nnz * sizeof(int));
    assert(cuda_error == cudaSuccess);

    cuda_error = cudaMemcpy(this->d_col, A.columns_indices().data(),
               this->nnz*sizeof(int), cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemcpy(this->d_row, A.row_extents().data(),
               (this->N+1)*sizeof(int), cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMemcpy(this->d_val, A.values().data(),
               this->nnz*sizeof(T), cudaMemcpyHostToDevice);
    assert(cuda_error == cudaSuccess);

    /******************************************/
    /* STEP 1: CREATE DESCRIPTORS FOR L AND U */
    /******************************************/
    cusparseCreateMatDescr(&descr_L);
    cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_UNIT);

    cusparseCreateMatDescr(&descr_U);
    cusparseSetMatType(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr_U, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatFillMode(descr_U, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatDiagType(descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);

    csrilu02Info_t info_A = 0;
    /**************************************************************************************************/
    /* STEP 2: QUERY HOW MUCH MEMORY USED IN LU FACTORIZATION AND THE TWO FOLLOWING SYSTEM INVERSIONS */
    /**************************************************************************************************/

    (cusparseCreateCsrilu02Info(&info_A));
    (cusparseCreateCsrsv2Info(&info_L));
    (cusparseCreateCsrsv2Info(&info_U));

    int pBufferSize_M, pBufferSize_L, pBufferSize_U;
    cusparseStatus = (cusparseDcsrilu02_bufferSize(this->cusparseHandle,
                                  this->N, this->nnz,
                                  this->description,
                                  this->d_val,
                                  this->d_row,
                                  this->d_col,
                                  info_A,
                                  &pBufferSize_M));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = (cusparseDcsrsv2_bufferSize(this->cusparseHandle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                this->N,
                                this->nnz, descr_L,
                                this->d_val,
                                this->d_row,
                                this->d_col, info_L,
                                &pBufferSize_L));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = (cusparseDcsrsv2_bufferSize(this->cusparseHandle,
                                CUSPARSE_OPERATION_NON_TRANSPOSE,
                                this->N,
                                this->nnz, descr_U,
                                this->d_val,
                                this->d_row,
                                this->d_col,
                                info_U,
                                &pBufferSize_U));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    int pBufferSize = max(pBufferSize_M, max(pBufferSize_L, pBufferSize_U));
    cuda_error = (cudaMalloc((void**)&pBuffer, pBufferSize));
    assert(cuda_error == cudaSuccess);

    /************************************************************************************************/
    /* STEP 3: ANALYZE THE THREE PROBLEMS: LU FACTORIZATION AND THE TWO FOLLOWING SYSTEM INVERSIONS */
    /************************************************************************************************/;
    cusparseStatus = (cusparseDcsrilu02_analysis(this->cusparseHandle,
                                this->N,
                                this->nnz,
                                this->description,
                                this->d_val,
                                this->d_row,
                                this->d_col,
                                info_A, CUSPARSE_SOLVE_POLICY_NO_LEVEL,
                                pBuffer));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    int structural_zero;
    cusparseStatus = cusparseXcsrilu02_zeroPivot(
            this->cusparseHandle, info_A, &structural_zero);
    //assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = (cusparseDcsrsv2_analysis(this->cusparseHandle,
                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                              this->N, this->nnz, descr_L,
                              this->d_val,
                              this->d_row, this->d_col,
                              info_L,
                              CUSPARSE_SOLVE_POLICY_NO_LEVEL, pBuffer));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseStatus = (cusparseDcsrsv2_analysis(this->cusparseHandle,
                              CUSPARSE_OPERATION_NON_TRANSPOSE,
                              this->N, this->nnz, descr_U,
                              this->d_val,
                              this->d_row,
                              this->d_col, info_U,
                              CUSPARSE_SOLVE_POLICY_USE_LEVEL, pBuffer));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /************************************/
    /* STEP 4: FACTORIZATION: A = L * U */
    /************************************/

    int numerical_zero;
    cusparseStatus = cusparseDcsrilu02(this->cusparseHandle,
                      this->N,
                      this->nnz,
                      this->description,
                      this->d_val,
                      this->d_row,
                      this->d_col, info_A,
                      CUSPARSE_SOLVE_POLICY_NO_LEVEL,
                      pBuffer);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);
    cusparseStatus = cusparseXcsrilu02_zeroPivot(
            this->cusparseHandle, info_A, &numerical_zero);
    //assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    (cudaMalloc(&d_z, this->N * sizeof(T)));
    (cudaMalloc(&d_x, this->N * sizeof(T)));
    (cudaMalloc(&d_y, this->N * sizeof(T)));
}

template<>
void GPUSparseLULinearSolver<float>::PreSolve(const SparseMatrixCSR<float>& A){
    throw std::invalid_argument(
            "GPUSparseLULinearSolver<float>::PreSolve not implemented");
}

template<class T>
void GPUSparseLULinearSolver<T>::Solve(const T* b, T* x){
    /*********************/
    /* STEP 5: L * z = x */
    /*********************/;

    (cudaMemcpy(d_x, b, this->N * sizeof(T),
                cudaMemcpyHostToDevice));

    cusparseStatus_t cusparseStatus;

    const T alpha = 1.0;
    cusparseStatus = (cusparseDcsrsv2_solve(this->cusparseHandle,
                           CUSPARSE_OPERATION_NON_TRANSPOSE,
                           this->N,
                           this->nnz,
                           &alpha, descr_L,
                           this->d_val,
                           this->d_row,
                           this->d_col,
                           info_L,
                           d_x,
                           d_z,
                           CUSPARSE_SOLVE_POLICY_NO_LEVEL,
                           pBuffer));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /*********************/
    /* STEP 5: U * y = z */
    /*********************/

    cusparseStatus = (cusparseDcsrsv2_solve(
            this->cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            this->N,
            this->nnz,
            &alpha,
            descr_U,
            this->d_val,
            this->d_row,
            this->d_col,
            info_U,
            d_z,
            d_y,
            CUSPARSE_SOLVE_POLICY_USE_LEVEL,
            pBuffer));
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    /********************************/
    /* MOVE THE RESULTS TO THE HOST */
    /********************************/
    cudaMemcpy(x, d_y, this->N * sizeof(T), cudaMemcpyDeviceToHost);
}

template<>
void GPUSparseLULinearSolver<float>::Solve(const float* b, float* x){
    throw std::invalid_argument(
            "GPUSparseLULinearSolver<float>::Solve not implemente");
}

template<class T>
void GPUSparseLULinearSolver<T>::Terminate(){
    if(this->pre_solved_){
        if(descr_L)
            cusparseDestroyMatDescr(descr_L);
        if(descr_U)
            cusparseDestroyMatDescr(descr_U);
        if(this->description)
            cusparseDestroyMatDescr(this->description);

        if(this->cusparseHandle)
            cusparseDestroy(this->cusparseHandle);
        if(this->d_col)
            cudaFree(this->d_col);
        if(this->d_row)
            cudaFree(this->d_row);
        if(this->d_val)
            cudaFree(this->d_val);

        if(d_z)
            cudaFree(d_z);
        if(d_x)
            cudaFree(d_x);
        if(d_y)
            cudaFree(d_y);
    }
}

template
class GPUSparseLULinearSolver<float>;
template
class GPUSparseLULinearSolver<double>;

}
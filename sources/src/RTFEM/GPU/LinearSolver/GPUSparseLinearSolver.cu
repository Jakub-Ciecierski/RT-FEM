#include "RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <assert.h>

namespace rtfem {

template<class T>
void GPUSparseLinearSolver<T>::Solve(){
    int N = 0, nz = 0, *I = NULL, *J = NULL;
    double *val = NULL;
    const double tol = 1e-5f;
    const int max_iter = 10000;
    double *x;
    double *rhs;
    double a, b, na, r0, r1;
    int *d_col, *d_row;
    double *d_val, *d_x, dot;
    double *d_r, *d_p, *d_Ax;
    int k;
    double alpha, beta, alpham1;

    /* Generate a random tridiagonal symmetric matrix in CSR format */
    N = 1048576;
    nz = (N-2)*3 + 4;
    I = (int *)malloc(sizeof(int)*(N+1));
    J = (int *)malloc(sizeof(int)*nz);
    val = (double *)malloc(sizeof(double)*nz);
    genTridiag(I, J, val, N, nz);

    x = (double *)malloc(sizeof(double)*N);
    rhs = (double *)malloc(sizeof(double)*N);

    for (int i = 0; i < N; i++)
    {
        rhs[i] = 1.0;
        x[i] = 0.0;
    }

    /* Get handle to the CUBLAS context */
    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;
    cudaError_t cuda_error;
    cublasStatus = cublasCreate(&cublasHandle);

    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    /* Get handle to the CUSPARSE context */
    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseHandle);
    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseMatDescr_t descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);

    assert(cusparseStatus == CUSPARSE_STATUS_SUCCESS);

    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    cuda_error = cudaMalloc((void **)&d_col, nz*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_row, (N+1)*sizeof(int));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_val, nz*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_x, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_r, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_p, N*sizeof(double));
    assert(cuda_error == cudaSuccess);
    cuda_error = cudaMalloc((void **)&d_Ax, N*sizeof(double));
    assert(cuda_error == cudaSuccess);

    cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_val, val, nz*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, rhs, N*sizeof(double), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz,
                   &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasDaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasDscal(cublasHandle, N, &b, d_p, 1);
            cublasStatus = cublasDaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasDcopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N,
                       N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublasStatus = cublasDdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasDaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasDaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasDdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaThreadSynchronize();
        printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }

    cudaMemcpy(x, d_x, N*sizeof(double), cudaMemcpyDeviceToHost);

    double rsum, diff, err = 0.0;

    for (int i = 0; i < N; i++)
    {
        rsum = 0.0;

        for (int j = I[i]; j < I[i+1]; j++)
        {
            rsum += val[j]*x[J[j]];
        }

        diff = fabs(rsum - rhs[i]);

        if (diff > err)
        {
            err = diff;
        }
    }

    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);

    free(I);
    free(J);
    free(val);
    free(x);
    free(rhs);
    cudaFree(d_col);
    cudaFree(d_row);
    cudaFree(d_val);
    cudaFree(d_x);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ax);

    printf("Test Summary:  Error amount = %f\n", err);
    return;
}

template<class T>
void GPUSparseLinearSolver<T>::genTridiag(int *I, int *J, double *val, int N, int nz)
{
    I[0] = 0, J[0] = 0, J[1] = 1;
    val[0] = (double)rand()/RAND_MAX + 10.0f;
    val[1] = (double)rand()/RAND_MAX;
    int start;

    for (int i = 1; i < N; i++)
    {
        if (i > 1)
        {
            I[i] = I[i-1]+3;
        }
        else
        {
            I[1] = 2;
        }

        start = (i-1)*3 + 2;
        J[start] = i - 1;
        J[start+1] = i;

        if (i < N-1)
        {
            J[start+2] = i + 1;
        }

        val[start] = val[start-1];
        val[start+1] = (double)rand()/RAND_MAX + 10.0f;

        if (i < N-1)
        {
            val[start+2] = (double)rand()/RAND_MAX;
        }
    }

    I[N] = nz;
}

template
class GPUSparseLinearSolver<float>;
template
class GPUSparseLinearSolver<double>;

}
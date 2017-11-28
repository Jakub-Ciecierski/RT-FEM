#include <RTFEM/GPU/VectorAdd/VectorAdd.cuh>

#include <stdio.h>
#include <cstdlib>

namespace rtfem {

#define SIZE 1024

CUDA_HOST_MEMBER
void VectorAddition::RunVectorAddTest() {
    int *a, *b, *c;

    auto size = SIZE * sizeof(int);
    a = (int *) malloc(size);
    b = (int *) malloc(size);
    c = (int *) malloc(size);
    for (int i = 0; i < SIZE; ++i) {
        a[i] = i;
        b[i] = i;
        c[i] = 0;
    }

    int *cuda_a, *cuda_b, *cuda_c;
    cudaError_t cuda_error;

    cuda_error = cudaMalloc(&cuda_a, size);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    cuda_error = cudaMalloc(&cuda_b, size);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    cuda_error = cudaMalloc(&cuda_c, size);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    cuda_error = cudaMemcpy(cuda_a, a, size, cudaMemcpyHostToDevice);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    cuda_error = cudaMemcpy(cuda_b, b, size, cudaMemcpyHostToDevice);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    cuda_error = cudaMemcpy(cuda_c, c, size, cudaMemcpyHostToDevice);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    VectorAdd << < 1, SIZE >> > (cuda_a, cuda_b, cuda_c, SIZE);
    cudaDeviceSynchronize();

    cuda_error = cudaMemcpy(c, cuda_c, size, cudaMemcpyDeviceToHost);
    if (cuda_error != cudaSuccess) {
        C_ERR(cuda_error);
    }

    for (int i = 0; i < 10; ++i)
        printf("c[%d] = %d\n", i, c[i]);

    free(a);
    free(b);
    free(c);

    cudaFree(cuda_a);
    cudaFree(cuda_a);
    cudaFree(cuda_a);
}

CUDA_KERNEL void VectorAdd(int *a, int *b, int *c, int n) {
    int i = threadIdx.x;
    if (i < n)
        c[i] = a[i] + b[i];
}

}
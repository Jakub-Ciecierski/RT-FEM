#ifndef _VECTOR_ADD_CUH_
#define _VECTOR_ADD_CUH_

#include <RTFEM/GPU/CudaFunctionMacros.cuh>

namespace rtfem {

class VectorAddition {
public:
    VectorAddition() = default;

    CUDA_HOST_MEMBER void RunVectorAddTest();
};

CUDA_KERNEL void VectorAdd(int *a, int *b, int *c, int n);

}

#endif
#ifndef PROJECT_CUDA_FUNCTION_MACROS_H
#define PROJECT_CUDA_FUNCTION_MACROS_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#define CUDA_HOST_MEMBER __host__
#define CUDA_KERNEL __global__
#else
#define CUDA_CALLABLE_MEMBER
#define CUDA_HOST_MEMBER
#define CUDA_KERNEL
#endif

#define C_ERR(err) (printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__),\
exit(EXIT_FAILURE))

#endif //PROJECT_CUDA_FUNCTION_MACROS_H

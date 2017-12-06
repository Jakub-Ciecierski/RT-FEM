#ifndef PROJECT_GPUMMMULTIPLICATION_H
#define PROJECT_GPUMMMULTIPLICATION_H

namespace rtfem {

enum class MatrixOperation{
    None, Transpose
};

template<class T>
class GPUMMMultiplication {
public:

    GPUMMMultiplication();

    ~GPUMMMultiplication();

    /**
     * A - m x k
     * B - k x n
     * C - m x n
     *
     * C = alpha * MatrixOperation(A)MatrixOperation(B) + beta*C
     */
    void Solve(const T* A, const T* B, T* C,
               T alpha, T beta,
               int m, int k, int n,
               MatrixOperation A_operation,
               MatrixOperation B_operation);
private:
};
}


#endif //PROJECT_GPUMMMULTIPLICATION_H

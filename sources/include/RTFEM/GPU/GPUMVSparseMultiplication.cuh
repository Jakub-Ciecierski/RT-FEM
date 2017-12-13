#ifndef PROJECT_GPUMVSPARSEMULTIPLICATION_H
#define PROJECT_GPUMVSPARSEMULTIPLICATION_H

namespace rtfem {

template<class T>
class SparseMatrixCSR;

template<class T>
class GPUMVSparseMultiplication {
public:

    GPUMVSparseMultiplication() = default;
    ~GPUMVSparseMultiplication() = default;

    void Solve(const SparseMatrixCSR<T>& A,
               T* x, T alpha,
               T* y, T beta);
private:
};
}


#endif //PROJECT_GPUMVSPARSEMULTIPLICATION_H

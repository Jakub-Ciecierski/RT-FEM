#ifndef PROJECT_DENSE2SPARSEMATRIX_H
#define PROJECT_DENSE2SPARSEMATRIX_H

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/SparseMatrixCSR.h>

namespace rtfem {


template<class T>
class Dense2SparseMatrix {
public:

    Dense2SparseMatrix() = default;
    ~Dense2SparseMatrix() = default;

    SparseMatrixCSR<T> Transform(
            const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&
            dense_matrix,
            MatrixType type);
private:
};
}


#endif //PROJECT_DENSE2SPARSEMATRIX_H

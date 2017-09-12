#ifndef PROJECT_MATRIXMATH_H
#define PROJECT_MATRIXMATH_H

#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/DataTypes.h>

namespace rtfem {

class MatrixMath {
public:
    MatrixMath();

    ~MatrixMath();

    /**
     * Computes determinant of a 2x2 matrix
     * @param matrix
     * @return
     */
    double ComputeDeterminant2(const Matrix &matrix);

    /**
     * Computes Determinant of Square matrix
     * @param matrix
     * @return
     */
    double ComputeDeterminant(const Matrix &matrix);

    /**
     * Removes specified row and column thus new matrix has dimension [N-1 x M-1]
     * @param matrix
     * @param row
     * @param column
     * @return
     */
    Matrix ContractMatrix(const Matrix &matrix,
                          unsigned int row,
                          unsigned int column);

    Matrix Transpose(const Matrix &matrix);
};
}

#endif //PROJECT_MATRIXMATH_H

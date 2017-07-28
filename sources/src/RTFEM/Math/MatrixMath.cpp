#include <stdexcept>
#include "RTFEM/Math/MatrixMath.h"

namespace rtfem {

MatrixMath::MatrixMath() {}

MatrixMath::~MatrixMath() {}

double MatrixMath::ComputeDeterminant2(const Matrix &matrix) {
    if(matrix.dimensions().row_count != 2 && matrix.dimensions().column_count != 2)
        throw std::invalid_argument("ComputeDeterminant2: Matrix must be 2x2");
    double a = matrix[0][0]* matrix[1][1];
    double b = matrix[0][1]* matrix[1][0];

    return a - b;
}

double MatrixMath::ComputeDeterminant(const Matrix &matrix) {
    unsigned int row_count =  matrix.dimensions().row_count;
    unsigned int column_count =  matrix.dimensions().column_count;
    if(column_count != row_count)
        throw std::invalid_argument("ComputeDeterminant: Requires square matrix");

    if(row_count == 2)
        return ComputeDeterminant2(matrix);

    int even_counter = 1;
    double determinant = 0;
    for(unsigned int i = 0; i < row_count; i++){
        auto a = matrix[0][i];
        auto contracted_matrix = ContractMatrix(matrix, 0, i);

        determinant += even_counter * a * ComputeDeterminant(contracted_matrix);

        even_counter = -even_counter;
    }

    return determinant;
}

Matrix MatrixMath::ContractMatrix(const Matrix& matrix, unsigned int row, unsigned int column){
    unsigned int row_count = matrix.dimensions().row_count;
    unsigned int column_count = matrix.dimensions().column_count;
    if(row >= row_count || column >= column_count)
        throw std::out_of_range("ContractMatrix: Index out of bounds");
    if(row_count == 1 || column_count == 1)
        throw std::invalid_argument("ContractMatrix: Matrix too small to contract");

    Matrix contracted_matrix(row_count-1, column_count-1);
    for(unsigned int i = 0; i < row_count; i++){
        if(i == row)
            continue;
        for(unsigned int j = 0; j < column_count; j++){
            if(j == column)
                continue;
            unsigned int contracted_i = i>row ? i-1 : i;
            unsigned int contracted_j = j>column? j-1 : j;
            contracted_matrix[contracted_i][contracted_j] = matrix[i][j];
        }
    }

    return contracted_matrix;
}

Matrix MatrixMath::Transpose(const Matrix &matrix) {
    Matrix m(matrix.dimensions().column_count,
             matrix.dimensions().row_count);

    for(unsigned int i = 0; i < matrix.dimensions().row_count; i++){
        for(unsigned int j = 0; j < matrix.dimensions().column_count; j++){
            m[j][i] = matrix[i][j];
        }
    }
    return m;
}

}
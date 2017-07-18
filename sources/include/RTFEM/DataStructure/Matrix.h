#ifndef PROJECT_MATRIX_H
#define PROJECT_MATRIX_H

#include <RTFEM/DataTypes.h>

#include <vector>
#include <iostream>

namespace rtfem {

struct MatrixDimension{
    MatrixDimension(UInt row_count, UInt column_count) :
            row_count(row_count), column_count(column_count){}

    UInt row_count;
    UInt column_count;
};

/**
 * Matrix data structure.
 */
class Matrix {
public:
    /**
     * Creates Matrix with row_count rows and column_count columns, initiated with zeros
     * @param row_count
     * @param column_count
     */
    Matrix(UInt row_count, UInt column_count);

    /**
     * Creates Matrix with row_count rows and column_count columns, initiated with value parameter
     * @param value
     * @param row_count
     * @param column_count
     */
    Matrix(Float value, UInt row_count, UInt column_count);

    Matrix(const MatrixDimension&& matrix_dimension);

    ~Matrix();

    const MatrixDimension& dimensions() const {return dimensions_;}

    /**
     * Returns i-th row.
     * @param i
     * @return
     */
    std::vector<Float>& operator[] (UInt i);
    const std::vector<Float>& operator[] (UInt i) const;

    /**
     * Returns i-th row.
     * @param i
     * @return
     */
    std::vector<Float>& GetRow(UInt i);
    const std::vector<Float>& GetRow(UInt i) const;

    bool operator==(const Matrix &rhs) const;
    bool operator!=(const Matrix &rhs) const;
private:
    void InitData(Float value);

    std::vector<std::vector<Float>> data_;

private:
    MatrixDimension dimensions_;
};

std::ostream& operator<<(std::ostream&, const Matrix&);

Matrix operator*(Float x, const Matrix& m);
Matrix operator*(const Matrix& m, Float x);

Matrix operator/(const Matrix& m, Float x);

Matrix operator*(const Matrix& m1, const Matrix& m2);
Matrix operator+(const Matrix& m1, const Matrix& m2);

}

#endif //PROJECT_MATRIX_H

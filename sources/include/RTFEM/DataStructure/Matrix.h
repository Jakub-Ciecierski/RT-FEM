#ifndef PROJECT_MATRIX_H
#define PROJECT_MATRIX_H

#include <RTFEM/DataTypes.h>

#include <vector>
#include <iostream>

namespace rtfem {

struct MatrixDimension {
    MatrixDimension(unsigned int row_count, unsigned int column_count) :
        row_count(row_count), column_count(column_count) {}

    unsigned int row_count;
    unsigned int column_count;
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
    Matrix(unsigned int row_count, unsigned int column_count);

    /**
     * Creates Matrix with row_count rows and column_count columns, initiated with value parameter
     * @param value
     * @param row_count
     * @param column_count
     */
    Matrix(double value, unsigned int row_count, unsigned int column_count);

    /**
     * Creates Matrix with elements taken from the initializer list.
     * The List contains rows (lists).
     *
     * Number of columns is the maximum count of all rows.
     * If any column has less elements then maximum, they are initialized with zero.
     *
     * @param lists
     */
    Matrix(std::initializer_list<std::initializer_list<double>> lists);

    Matrix(const MatrixDimension &&matrix_dimension);

    ~Matrix();

    const MatrixDimension &dimensions() const { return dimensions_; }

    /**
     * Returns i-th row.
     * @param i
     * @return
     */
    std::vector<double> &operator[](unsigned int i);
    const std::vector<double> &operator[](unsigned int i) const;

    /**
     * Returns i-th row.
     * @param i
     * @return
     */
    std::vector<double> &GetRow(unsigned int i);
    const std::vector<double> &GetRow(unsigned int i) const;

    bool operator==(const Matrix &rhs) const;
    bool operator!=(const Matrix &rhs) const;

    Matrix &operator+=(const Matrix &rhs);

private:
    void InitData(double value);

    std::vector<std::vector<double>> data_;

private:
    MatrixDimension dimensions_;
};

std::ostream &operator<<(std::ostream &, const Matrix &);

Matrix operator*(double x, const Matrix &m);
Matrix operator*(const Matrix &m, double x);

Matrix operator/(const Matrix &m, double x);

Matrix operator*(const Matrix &m1, const Matrix &m2);
Matrix operator+(const Matrix &m1, const Matrix &m2);

}

#endif //PROJECT_MATRIX_H

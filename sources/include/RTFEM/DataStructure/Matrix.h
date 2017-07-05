#ifndef PROJECT_MATRIX_H
#define PROJECT_MATRIX_H

#include <RTFEM/DataTypes.h>

#include <vector>

namespace rtfem {

struct MatrixDimension{
    MatrixDimension(UInt row_count, UInt column_count) :
            row_count(row_count), column_count(column_count){}

    UInt row_count;
    UInt column_count;
};

class Matrix {
public:
    /**
     * Creates Matrix with row_count rows and column_count columns
     * @param row_count
     * @param column_count
     */
    Matrix(UInt row_count, UInt column_count);
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

private:
    void InitData();

    std::vector<std::vector<Float>> data_;

    MatrixDimension dimensions_;
};

}

#endif //PROJECT_MATRIX_H

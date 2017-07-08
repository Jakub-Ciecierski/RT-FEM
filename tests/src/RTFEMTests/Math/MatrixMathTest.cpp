#include "RTFEMTests/Math/MatrixMathTest.h"

#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/Math/MatrixMath.h>
#include <RTFEM/DataTypes.h>

void MatrixMathTest::SetUp() {
    matrix_ = rtfem::make_unique<rtfem::Matrix>(3,3);

    (*matrix_)[0][0] = 0;
    (*matrix_)[0][1] = 1;
    (*matrix_)[0][2] = 2;

    (*matrix_)[1][0] = 3;
    (*matrix_)[1][1] = 4;
    (*matrix_)[1][2] = 5;

    (*matrix_)[2][0] = 6;
    (*matrix_)[2][1] = 7;
    (*matrix_)[2][2] = 8;

    matrix_math_ = rtfem::make_unique<rtfem::MatrixMath>();
}

void MatrixMathTest::TearDown() {
}

TEST_F(MatrixMathTest, ContractMatrix_RowIndexOutOfBounds_ThrownOutOfRangeExpection){
    EXPECT_THROW(matrix_math_->ContractMatrix(std::move(*matrix_), 4, 2), std::out_of_range);
}

TEST_F(MatrixMathTest, ContractMatrix_ColumnIndexOutOfBounds_ThrownOutOfRangeExpection){
    EXPECT_THROW(matrix_math_->ContractMatrix(std::move(*matrix_), 2, 4), std::out_of_range);
}

TEST_F(MatrixMathTest, ContractMatrix_MatrixWithSingleRow_ThrownInvalidArgumentExpection){
    auto small_matrix = rtfem::make_unique<rtfem::Matrix>(1,1);
    EXPECT_THROW(matrix_math_->ContractMatrix(std::move(*small_matrix), 0, 0), std::invalid_argument);
}

TEST_F(MatrixMathTest, ContractMatrix_Matrix_DimensionsSmallerByOne){
    auto matrix = matrix_math_->ContractMatrix(std::move(*matrix_), 1, 1);

    const rtfem::UInt expected_row_count = 2;
    const rtfem::UInt expected_column_count = 2;

    EXPECT_EQ(expected_row_count, matrix.dimensions().row_count);
    EXPECT_EQ(expected_column_count, matrix.dimensions().column_count);
}

TEST_F(MatrixMathTest, ContractMatrix_Matrix_ProperValues){
    auto contracted_matrix = matrix_math_->ContractMatrix(std::move(*matrix_), 1, 1);

    EXPECT_EQ((rtfem::Float)0.0, contracted_matrix[0][0]);
    EXPECT_EQ((rtfem::Float)2.0, contracted_matrix[0][1]);

    EXPECT_EQ((rtfem::Float)6.0, contracted_matrix[1][0]);
    EXPECT_EQ((rtfem::Float)8.0, contracted_matrix[1][1]);
}

TEST_F(MatrixMathTest, ComputeDeterminant2_ProperValue){
    rtfem::Matrix matrix(2,2);
    matrix[0][0] = 4;
    matrix[0][1] = 6;
    matrix[1][0] = 3;
    matrix[1][1] = 8;

    const rtfem::Float expected_det = 14;
    rtfem::Float det = matrix_math_->ComputeDeterminant2(matrix);

    EXPECT_EQ(expected_det, det);
}

TEST_F(MatrixMathTest, ComputeDeterminant_Matrix3ProperValue){
    rtfem::Matrix matrix(3,3);
    matrix[0][0] = 6;
    matrix[0][1] = 1;
    matrix[0][2] = 1;

    matrix[1][0] = 4;
    matrix[1][1] = -2;
    matrix[1][2] = 5;

    matrix[2][0] = 2;
    matrix[2][1] = 8;
    matrix[2][2] = 7;

    const rtfem::Float expected_det = -306;
    rtfem::Float det = matrix_math_->ComputeDeterminant(matrix);

    EXPECT_EQ(expected_det, det);
}

TEST_F(MatrixMathTest, ComputeDeterminant_Matrix4ProperValue){
    rtfem::Matrix matrix(4,4);
    matrix[0][0] = 11;
    matrix[0][1] = 223;
    matrix[0][2] = 4;
    matrix[0][3] = 54;

    matrix[1][0] = 5;
    matrix[1][1] = 5;
    matrix[1][2] = 58;
    matrix[1][3] = 6;

    matrix[2][0] = 8;
    matrix[2][1] = 8;
    matrix[2][2] = 7;
    matrix[2][3] = 8;

    matrix[3][0] = 8;
    matrix[3][1] = 8;
    matrix[3][2] = 8;
    matrix[3][3] = 9;

    const rtfem::Float expected_det = 89252;
    rtfem::Float det = matrix_math_->ComputeDeterminant(matrix);

    EXPECT_EQ(expected_det, det);
}

TEST_F(MatrixMathTest, Transpose_ProperMatrix){
    rtfem::Matrix m(3, 2);
    m[0][0] = 1;
    m[0][1] = 2;

    m[1][0] = 3;
    m[1][1] = 4;

    m[2][0] = 5;
    m[2][1] = 6;

    rtfem::Matrix expected_m(2,3);
    expected_m[0][0] = 1;
    expected_m[0][1] = 3;
    expected_m[0][2] = 5;

    expected_m[1][0] = 2;
    expected_m[1][1] = 4;
    expected_m[1][2] = 6;

    EXPECT_EQ(expected_m, matrix_math_->Transpose(m));
}

TEST_F(MatrixMathTest, Transpose_TransposedTwice_BackToOriginal){
    rtfem::Matrix m(3, 2);
    m[0][0] = 1;
    m[0][1] = 2;

    m[1][0] = 3;
    m[1][1] = 4;

    m[2][0] = 5;
    m[2][1] = 6;

    EXPECT_EQ(m, matrix_math_->Transpose(matrix_math_->Transpose(m)));
}
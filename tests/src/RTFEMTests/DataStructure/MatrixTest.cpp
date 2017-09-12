#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/Memory/UniquePointer.h>

#include <RTFEMTests/DataStructure/MatrixTest.h>

void MatrixTest::SetUp() {
    small_matrix_.row_count = 5;
    small_matrix_.column_count = 3;
    small_matrix_.matrix =
        CreateMatrix(small_matrix_.row_count, small_matrix_.column_count);
}

void MatrixTest::TearDown() {

}

std::unique_ptr<rtfem::Matrix> MatrixTest::CreateMatrix(unsigned int r,
                                                        unsigned int c) {
    return rtfem::make_unique<rtfem::Matrix>(r, c);
}

// -------------------------------------

TEST_F(MatrixTest, Matrix_GetDimensions_ProperDimensions) {
    EXPECT_EQ(small_matrix().row_count,
              small_matrix().matrix->dimensions().row_count);
    EXPECT_EQ(small_matrix().column_count,
              small_matrix().matrix->dimensions().column_count);
}

TEST_F(MatrixTest, IndexOutOfBounds_ThrowsException) {
    auto i = small_matrix().row_count;

    EXPECT_THROW((*small_matrix().matrix)[i], std::out_of_range);
}

TEST_F(MatrixTest, IndexNotOutOfBounds_DoesNotThrowException) {
    (*small_matrix().matrix)[0];
}

TEST_F(MatrixTest, MatrixConstructed_FullOfZeros) {
    for (unsigned int i = 0; i > small_matrix().matrix->dimensions().row_count;
         i++) {
        for (unsigned int j = 0;
             j > small_matrix().matrix->dimensions().column_count; j++) {
            EXPECT_EQ(0, (*small_matrix().matrix)[i][j]);
        }
    }
}

TEST_F(MatrixTest, SetValue_GetSameValue) {
    const double x = 5.0;

    (*small_matrix().matrix)[0][0] = x;

    EXPECT_EQ(x, (*small_matrix().matrix)[0][0]);
}

TEST_F(MatrixTest, MatrixEqualityOperator_TwoEqualMatrices_ShouldBeEqual) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(3, 2);
    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[1][0] = 3;
    m2[1][1] = 4;
    m2[2][0] = 5;
    m2[2][1] = 6;

    EXPECT_EQ(true, m1 == m2);
}

TEST_F(MatrixTest,
       MatrixEqualityOperator_TwoNotEqualMatricesDifferentDimension_ShouldNotBeEqual) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(2, 4);

    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[0][2] = 3;
    m2[0][3] = 4;
    m2[1][0] = 5;
    m2[1][1] = 6;
    m2[1][2] = 7;
    m2[1][3] = 8;

    EXPECT_EQ(false, m1 == m2);
}

TEST_F(MatrixTest,
       MatrixEqualityOperator_TwoNotEqualMatricesSameDimension_ShouldNotBeEqual) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(3, 2);
    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[1][0] = 0;
    m2[1][1] = 4;
    m2[2][0] = 5;
    m2[2][1] = 6;

    EXPECT_EQ(false, m1 == m2);
}

TEST_F(MatrixTest, MatrixScalarMultiplication_RightSide) {
    const double x = 5.0;
    const double expected_value = 25;

    (*small_matrix().matrix)[0][0] = x;

    auto multiplied_matrix = (*small_matrix().matrix) * x;

    EXPECT_EQ(expected_value, multiplied_matrix[0][0]);
}

TEST_F(MatrixTest, MatrixScalarMultiplication_LeftSide) {
    const double x = 5.0;
    const double expected_value = 25;

    (*small_matrix().matrix)[0][0] = x;

    auto multiplied_matrix = x * (*small_matrix().matrix);

    EXPECT_EQ(expected_value, multiplied_matrix[0][0]);
}

TEST_F(MatrixTest, MatrixScalarDivision_RightSide) {
    const double expected_value = 5;

    (*small_matrix().matrix)[0][0] = 10.0;

    auto multiplied_matrix = (*small_matrix().matrix) / 2.0;

    EXPECT_EQ(expected_value, multiplied_matrix[0][0]);
}

TEST_F(MatrixTest,
       MatrixMatrixMultiplication_IncorrectDimensions_ExpectionThrown) {
    rtfem::Matrix m1(2, 3);
    rtfem::Matrix m2(4, 5);
    EXPECT_THROW(m1 * m2, std::invalid_argument);
}

TEST_F(MatrixTest, MatrixMatrixMultiplication_CorrectAnswer) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(2, 4);

    m2[0][0] = 1;
    m2[0][1] = 2;
    m2[0][2] = 3;
    m2[0][3] = 4;
    m2[1][0] = 5;
    m2[1][1] = 6;
    m2[1][2] = 7;
    m2[1][3] = 8;

    auto m = m1 * m2;

    rtfem::Matrix expected_result
        (m1.dimensions().row_count, m2.dimensions().column_count);
    expected_result[0][0] = 11;
    expected_result[0][1] = 14;
    expected_result[0][2] = 17;
    expected_result[0][3] = 20;

    expected_result[1][0] = 23;
    expected_result[1][1] = 30;
    expected_result[1][2] = 37;
    expected_result[1][3] = 44;

    expected_result[2][0] = 35;
    expected_result[2][1] = 46;
    expected_result[2][2] = 57;
    expected_result[2][3] = 68;

    EXPECT_EQ(expected_result, m);
}

TEST_F(MatrixTest, MatrixAddition_IncorrectDimensions_ExceptionThrown) {
    rtfem::Matrix m1(3, 2);
    rtfem::Matrix m2(2, 4);

    EXPECT_THROW(m1 + m2, std::invalid_argument);
}

TEST_F(MatrixTest, MatrixAddition_ProperResult) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(3, 2);
    m2[0][0] = 10;
    m2[0][1] = 2;
    m2[1][0] = 4;
    m2[1][1] = 7;
    m2[2][0] = 8;
    m2[2][1] = 1;

    rtfem::Matrix expected_result(3, 2);
    expected_result[0][0] = 11;
    expected_result[0][1] = 4;
    expected_result[1][0] = 7;
    expected_result[1][1] = 11;
    expected_result[2][0] = 13;
    expected_result[2][1] = 7;

    EXPECT_EQ(expected_result, m1 + m2);
}

TEST_F(MatrixTest, MatrixReferenceAddition_ProperResult) {
    rtfem::Matrix m1(3, 2);
    m1[0][0] = 1;
    m1[0][1] = 2;
    m1[1][0] = 3;
    m1[1][1] = 4;
    m1[2][0] = 5;
    m1[2][1] = 6;

    rtfem::Matrix m2(3, 2);
    m2[0][0] = 10;
    m2[0][1] = 2;
    m2[1][0] = 4;
    m2[1][1] = 7;
    m2[2][0] = 8;
    m2[2][1] = 1;

    rtfem::Matrix expected_result(3, 2);
    expected_result[0][0] = 11;
    expected_result[0][1] = 4;
    expected_result[1][0] = 7;
    expected_result[1][1] = 11;
    expected_result[2][0] = 13;
    expected_result[2][1] = 7;

    m1 += m2;

    EXPECT_EQ(expected_result, m1);
}

TEST_F(MatrixTest, MatrixInitializedList_AllColumnsProvided_CorrectDimensions) {
    rtfem::Matrix m{{1, 2, 3}, {4, 5, 6}};
    constexpr unsigned int expected_row_count = 2;
    constexpr unsigned int expected_column_count = 3;

    EXPECT_EQ(expected_row_count, m.dimensions().row_count);
    EXPECT_EQ(expected_column_count, m.dimensions().column_count);
}

TEST_F(MatrixTest,
       MatrixInitializedList_NotAllColumnsProvided_CorrectDimensions) {
    rtfem::Matrix m{{1, 2, 3}, {4, 6}};
    unsigned int expected_row_count = 2;
    unsigned int expected_column_count = 3;

    EXPECT_EQ(expected_row_count, m.dimensions().row_count);
    EXPECT_EQ(expected_column_count, m.dimensions().column_count);
}

TEST_F(MatrixTest,
       MatrixInitializedList_NotAllColumnsProvided_InitializedWithZero) {
    rtfem::Matrix m{{1, 2, 3}, {4, 6}};

    EXPECT_EQ(0, m[1][2]);
}
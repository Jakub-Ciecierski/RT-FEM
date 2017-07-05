#include "gtest/gtest.h"

#include <RTFEM/DataTypes.h>
#include <RTFEM/DataStructure/Matrix.h>
#include <RTFEM/Memory/UniquePointer.h>

#include <RTFEMTests/DataStructure/MatrixTest.h>

void MatrixTest::SetUp(){
    small_matrix_.row_count = 5;
    small_matrix_.column_count = 3;
    small_matrix_.matrix = CreateMatrix(small_matrix_.row_count, small_matrix_.column_count);
}

void MatrixTest::TearDown(){

}

std::unique_ptr<rtfem::Matrix> MatrixTest::CreateMatrix(rtfem::UInt r, rtfem::UInt c){
    return rtfem::make_unique<rtfem::Matrix>(r,c);
}

// -------------------------------------

TEST_F(MatrixTest, Matrix_GetDimensions_ProperDimensions) {
    EXPECT_EQ(small_matrix().row_count, small_matrix().matrix->dimensions().row_count);
    EXPECT_EQ(small_matrix().column_count, small_matrix().matrix->dimensions().column_count);
}

TEST_F(MatrixTest, IndexOutOfBounds_ThrowsException) {
    auto i = small_matrix().row_count;

    EXPECT_THROW((*small_matrix().matrix)[i], std::out_of_range);
}

TEST_F(MatrixTest, IndexNotOutOfBounds_DoesNotThrowException) {
    (*small_matrix().matrix)[0];
}

TEST_F(MatrixTest, SetValue_GetSameValue) {
    const rtfem::Float x = 5.0;

    (*small_matrix().matrix)[0][0] = x;

    EXPECT_EQ(x, (*small_matrix().matrix)[0][0]);
}
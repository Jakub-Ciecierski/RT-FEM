#include "RTFEMTests/GPU/GPUMatrixMultiplicationTest.h"

#include "RTFEM/GPU/GPUMatrixMultiplication.cuh"
#include "RTFEM/DataTypes.h"

void GPUMatrixMultiplicationTest::SetUp() {
}

void GPUMatrixMultiplicationTest::TearDown() {
}

TEST_F(GPUMatrixMultiplicationTest, MatrixVector_ProperMultiplicationResults){
    constexpr int n = 3;
    Eigen::Matrix<double, n, n> A;
    A << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;

    Eigen::Vector<double, n> x;
    x << 1, 2, 3;

    Eigen::Vector<double, n> y;

    rtfem::GPUMatrixMultiplication<double> gpu_multiplication;
    gpu_multiplication.PreSolve(A.data(), n);
    gpu_multiplication.Solve(x.data(), 1, y.data(), 0);

    EXPECT_EQ(14, y(0));
    EXPECT_EQ(32, y(1));
    EXPECT_EQ(50, y(2));
}
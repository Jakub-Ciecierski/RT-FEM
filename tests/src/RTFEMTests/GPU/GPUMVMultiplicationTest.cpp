#include "RTFEMTests/GPU/GPUMVMultiplicationTest.h"

#include "RTFEM/GPU/GPUMVMultiplication.cuh"
#include "RTFEM/DataTypes.h"

void GPUMVMultiplicationTest::SetUp() {
}

void GPUMVMultiplicationTest::TearDown() {
}

TEST_F(GPUMVMultiplicationTest, MatrixVector_ProperMultiplicationResults){
    constexpr int n = 3;
    Eigen::Matrix<double, n, n> A;
    A << 1, 2, 3,
            4, 5, 6,
            7, 8, 9;

    Eigen::Vector<double, n> x;
    x << 1, 2, 3;

    Eigen::Vector<double, n> y;

    rtfem::GPUMVMultiplication<double> gpu_multiplication;
    gpu_multiplication.PreSolve(A.data(), n);
    gpu_multiplication.Solve(x.data(), 1, y.data(), 0);

    EXPECT_EQ(14, y(0));
    EXPECT_EQ(32, y(1));
    EXPECT_EQ(50, y(2));
}
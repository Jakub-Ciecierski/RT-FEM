#include "RTFEMTests/GPU/GPUMMMultiplicationTest.h"

#include "RTFEM/GPU/GPUMMMultiplication.cuh"
#include <RTFEM/DataTypes.h>

void GPUMMMultiplicationTest::SetUp() {
}

void GPUMMMultiplicationTest::TearDown() {
}

TEST_F(GPUMMMultiplicationTest, Multiplying2Matrices_ProperResult){
    constexpr int m = 3;
    constexpr int k = 2;
    constexpr int n = 4;

    Eigen::Matrix<double, m, k> A;
    A << 1, 2,
            3, 4,
            5, 6;
    Eigen::Matrix<double, k, n> B;
    B << 1, 2, 3, 4,
            5, 6, 7, 8;
    Eigen::Matrix<double, m, n> C;

    rtfem::GPUMMMultiplication<double> gpu_mm;
    gpu_mm.Solve(A.data(), B.data(), C.data(),
                 1, 0,
                 m, k, n,
                 rtfem::MatrixOperation::None,
                 rtfem::MatrixOperation::None);

    Eigen::Matrix<double, m, n> C_expected;
    C_expected << 11, 14, 17, 20,
            23, 30, 37, 44,
            35, 46, 57, 68;

    EXPECT_EQ(C_expected, C);
}

TEST_F(GPUMMMultiplicationTest, Multiplying2AndAddingMatrices_ProperResult){
    constexpr int m = 3;
    constexpr int k = 2;
    constexpr int n = 4;

    Eigen::Matrix<double, m, k> A;
    A << 1, 2,
            3, 4,
            5, 6;
    Eigen::Matrix<double, k, n> B;
    B << 1, 2, 3, 4,
            5, 6, 7, 8;
    Eigen::Matrix<double, m, n> C;
    C << 1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12;

    rtfem::GPUMMMultiplication<double> gpu_mm;
    gpu_mm.Solve(A.data(), B.data(), C.data(),
                 1, 1,
                 m, k, n,
                 rtfem::MatrixOperation::None,
                 rtfem::MatrixOperation::None);

    Eigen::Matrix<double, m, n> C_expected;
    C_expected << 12, 16, 20, 24,
            28, 36, 44, 52,
            44, 56, 68, 80;

    EXPECT_EQ(C_expected, C);
}

TEST_F(GPUMMMultiplicationTest, ASD){
    constexpr int m = 12;
    constexpr int k = 15;

    Eigen::Matrix<double, m, k> A;
    A <<    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 60, 60, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 7, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 20, 0, 6, 0, 0, 8, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0,
            7, 0, 0, 0, 0, 0, 0, 0, 60, 20, 0, 5, 0, 5, 0,
            0, 0, 2, 5, 6, 7, 0, 0, 0, 0, 0, 0, 2, 6, 0,
            6, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0,
            0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor> A_T = A.transpose();

    Eigen::Matrix<double, m, m> B;
    B <<    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 2, 0, 0, 7, 0, 0, 0, 0, 7, 0, 0,
            0, 0, 2, 8, 5, 4, 2, 0, 0, 8, 6, 0,
            0, 0, 9, 2, 0, 0, 6, 0, 0, 0, 0, 2,
            0, 0, 9, 0, 2, 0, 0, 0, 0, 0, 0, 0,
            0, 6, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
            0, 0, 5, 0, 2, 0, 2, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 2, 0, 4, 5, 0,
            0, 0, 3, 0, 2, 0, 0, 0, 2, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 6, 0, 0, 2, 0, 0,
            0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2;

    auto C_expected =  A_T * B * A;

    Eigen::Matrix<double, k, m> C1;
    rtfem::GPUMMMultiplication<double> gpu_mm;
    gpu_mm.Solve(A_T.data(), B.data(), C1.data(),
                 1, 0,
                 k, m, m,
                 rtfem::MatrixOperation::None,
                 rtfem::MatrixOperation::None);

    Eigen::Matrix<double, k, k> C2;
    gpu_mm.Solve(C1.data(), A.data(), C2.data(),
                 1, 0,
                 k, m, k,
                 rtfem::MatrixOperation::None,
                 rtfem::MatrixOperation::None);
    EXPECT_EQ(C_expected, C2);
}

TEST_F(GPUMMMultiplicationTest, ASD2){
    constexpr int m = 3;
    constexpr int k = 4;

    Eigen::Matrix<double, m, k> A;
    A <<    1, 2, 3, 4,
            5, 6, 7, 9,
            0, 0, 1, 0;

    Eigen::Matrix<double, k, m> A_T2;
    A_T2 <<  1, 5, 0,
            2, 6, 0,
            3, 7, 1,
            4, 9, 0;

    Eigen::Matrix<double, k, m, Eigen::ColMajor> A_T = A.transpose();

    Eigen::Matrix<double, m, m> B;
    B <<    2, 0, 0,
            0, 2, 0,
            0, 0, 2;

    auto C_expected1 = A_T * B;

    Eigen::Matrix<double, k, m> C1;
    rtfem::GPUMMMultiplication<double> gpu_mm;
    gpu_mm.Solve(A_T.data(), B.data(), C1.data(),
                 1, 0,
                 k, m, m,
                 rtfem::MatrixOperation::None,
                 rtfem::MatrixOperation::None);

    EXPECT_EQ(C_expected1, C1);

}

#include "RTFEMTests/GPU/GPUMVSparseMultiplicationTest.h"

#include <RTFEM/DataTypes.h>
#include "RTFEM/DataStructure/Dense2SparseMatrix.h"
#include "RTFEM/GPU/GPUMVSparseMultiplication.cuh"

void GPUMVSparseMultiplicationTest::SetUp() {
}

void GPUMVSparseMultiplicationTest::TearDown() {
}

TEST_F(GPUMVSparseMultiplicationTest, DenseMatrix){
    Eigen::Matrix<double, 4, 4> dense_matrix;
    dense_matrix <<
            0, 0, 0, 0,
            5, 8, 0, 0,
            0, 0, 3, 0,
            0, 6, 0, 0;
    Eigen::Vector<double, 4> vector;
    vector << 1, 2, 3, 4;

    rtfem::Dense2SparseMatrix<double> dense2sparse;
    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::General);

    Eigen::Vector<double, 4> out_vector;
    rtfem::GPUMVSparseMultiplication<double> gpu_mv_sparse;

    gpu_mv_sparse.PreSolve(sparse_matrix);
    gpu_mv_sparse.Solve(vector.data(), 1,
                        out_vector.data(), 0);

    Eigen::Vector<double, 4> expected_result;
    expected_result << 0, 21, 9, 12;

    EXPECT_EQ(expected_result, out_vector);
}

TEST_F(GPUMVSparseMultiplicationTest, SparseMatrix){
    Eigen::Matrix<double, 5, 5> dense_matrix;
    dense_matrix <<
             1,  -1,  0,  -3, 0,
            -1,  5,  0,  0,  0,
            0,   0,  4,  6,  4,
            -3,  0,  6,  7,  0,
            0,   0,  4,  0,  -5;
    Eigen::Vector<double, 5> vector;
    vector << 1, 2, 3, 4, 5;

    rtfem::Dense2SparseMatrix<double> dense2sparse;
    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::General);

    Eigen::Vector<double, 5> out_vector;
    rtfem::GPUMVSparseMultiplication<double> gpu_mv_sparse;

    gpu_mv_sparse.PreSolve(sparse_matrix);
    gpu_mv_sparse.Solve(vector.data(), 1,
                        out_vector.data(), 0);

    auto expected_result = dense_matrix * vector;

    std::cout << "expected_result" << std::endl;
    std::cout << expected_result << std::endl;
    std::cout << "out_vector" << std::endl;
    std::cout << out_vector << std::endl;

    EXPECT_EQ(expected_result, out_vector);
}
#include "RTFEMTests/DataStructure/Dense2SparseMatrixTest.h"

#include "RTFEM/DataStructure/Dense2SparseMatrix.h"

void Dense2SparseMatrixTest::SetUp() {
}

void Dense2SparseMatrixTest::TearDown() {
}

TEST_F(Dense2SparseMatrixTest, Transform_GeneralMatrix_1){
    rtfem::Dense2SparseMatrix<double> dense2sparse;

    Eigen::Matrix<double, 4, 4> dense_matrix;
    dense_matrix <<
                 0, 0, 0, 0,
            5, 8, 0, 0,
            0, 0, 3, 0,
            0, 6, 0, 0;

    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::General);

    std::vector<double> values = {5, 8, 3, 6};
    std::vector<int> row_extends = {0, 0, 2, 3, 4};
    std::vector<int> column_indices = {0, 1, 2, 1};

    EXPECT_EQ(values, sparse_matrix.values());
    EXPECT_EQ(row_extends, sparse_matrix.row_extents());
    EXPECT_EQ(column_indices, sparse_matrix.columns_indices());
}

TEST_F(Dense2SparseMatrixTest, Transform_GeneralMatrix_2){
    rtfem::Dense2SparseMatrix<double> dense2sparse;

    Eigen::Matrix<double, 4, 6> dense_matrix;
    dense_matrix <<
                 10, 20, 0, 0, 0, 0,
            0, 30, 0, 40, 0, 0,
            0, 0, 50, 60, 70, 0,
            0, 0, 0, 0, 0, 80;

    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::General);

    std::vector<double> values = {10, 20, 30, 40, 50, 60, 70, 80};
    std::vector<int> row_extends = {0, 2, 4, 7, 8};
    std::vector<int> column_indices = {0, 1, 1, 3, 2, 3, 4, 5};

    EXPECT_EQ(values, sparse_matrix.values());
    EXPECT_EQ(row_extends, sparse_matrix.row_extents());
    EXPECT_EQ(column_indices, sparse_matrix.columns_indices());
}

TEST_F(Dense2SparseMatrixTest, Transform_SymmetricMatrix){
    rtfem::Dense2SparseMatrix<double> dense2sparse;

    Eigen::Matrix<double, 5, 5> dense_matrix;
    dense_matrix <<
            1,  -1,  0,  -3, 0,
            -1,  5,  0,  0,  0,
            0,   0,  4,  6,  4,
            -3,  0,  6,  7,  0,
            0,   0,  4,  0,  -5;

    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::Symmetric);

    std::vector<double> values = {1, -1, -3, 5, 4, 6, 4, 7, -5};

    std::vector<int> row_extends = {0, 3, 4, 7, 8, 9};
    std::vector<int> column_indices = {0, 1, 3, 1, 2, 3, 4, 3, 4};

    EXPECT_EQ(values, sparse_matrix.values());
    EXPECT_EQ(row_extends, sparse_matrix.row_extents());
    EXPECT_EQ(column_indices, sparse_matrix.columns_indices());
}
#include "RTFEMTests/GPU/LinearSolver/GPUSparseLinearSolverTest.h"

#include <RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh>
#include <RTFEM/DataTypes.h>
#include "RTFEM/DataStructure/Dense2SparseMatrix.h"

void GPUSparseLinearSolverTest::SetUp() {
}

void GPUSparseLinearSolverTest::TearDown() {
}

TEST_F(GPUSparseLinearSolverTest, Integration){
    Eigen::Matrix<double, 5, 5> dense_matrix;
    dense_matrix <<
                 1,  -1,  0,  -3, 0,
            -1,  5,  0,  0,  0,
            0,   0,  4,  6,  4,
            -3,  0,  6,  7,  0,
            0,   0,  4,  0,  -5;
    Eigen::Vector<double, 5> b;
    b << -13, 9, 56, 43, -13;
    Eigen::Vector<double, 5> x;
    x << 0, 0, 0, 0, 0;

    rtfem::Dense2SparseMatrix<double> dense2sparse;
    auto sparse_matrix = dense2sparse.Transform(dense_matrix,
                                                rtfem::MatrixType::General);

    rtfem::GPUSparseLinearSolver<double> gpu_sparse_solver;
    gpu_sparse_solver.PreSolve(sparse_matrix);
    gpu_sparse_solver.Solve(b.data(), x.data());

}
#include "RTFEMTests/GPU/LinearSolver/GPUSparseLinearSolverTest.h"

#include <RTFEM/GPU/LinearSolver/GPUSparseCGLinearSolver.cuh>
#include <RTFEM/GPU/LinearSolver/GPUSparseCGPrecondLinearSolver.cuh>
#include <RTFEM/GPU/LinearSolver/GPUSparseLULinearSolver.cuh>

#include <RTFEM/DataTypes.h>
#include <RTFEM/Memory/UniquePointer.h>
#include "RTFEM/DataStructure/Dense2SparseMatrix.h"

void GPUSparseLinearSolverTest::SetUp() {
    dense_matrix <<
                 1,  -1,  0,  -3, 0,
            -1,  5,  0,  0,  0,
            0,   0,  4,  6,  4,
            -3,  0,  6,  7,  0,
            0,   0,  4,  0,  -5;
    b << -13, 9, 56, 43, -13;
    x << 0, 0, 0, 0, 0;

    rtfem::Dense2SparseMatrix<double> dense2sparse;
    auto sparse_matrix_local = dense2sparse.Transform(
            dense_matrix, rtfem::MatrixType::General);

    sparse_matrix = rtfem::make_unique<rtfem::SparseMatrixCSR<double>>(
            sparse_matrix_local.values(),
            sparse_matrix_local.row_extents(),
            sparse_matrix_local.columns_indices(),
            sparse_matrix_local.m(),
            sparse_matrix_local.n(),
            sparse_matrix_local.type()
    );
}

void GPUSparseLinearSolverTest::TearDown() {
}

TEST_F(GPUSparseLinearSolverTest, GPUSparseCGLinearSolver_Integration){
    rtfem::GPUSparseCGLinearSolver<double> gpu_sparse_solver;
    gpu_sparse_solver.PreSolve(*sparse_matrix);
    gpu_sparse_solver.Solve(b.data(), x.data());

    std::cout << "x" << std::endl;
    std::cout << x << std::endl;
}

TEST_F(GPUSparseLinearSolverTest, GPUSparseCGPrecondLinearSolver_Integration){
    rtfem::GPUSparseCGPrecondLinearSolver<double> gpu_sparse_solver;
    gpu_sparse_solver.PreSolve(*sparse_matrix);
    gpu_sparse_solver.Solve(b.data(), x.data());

    std::cout << "x" << std::endl;
    std::cout << x << std::endl;
}

TEST_F(GPUSparseLinearSolverTest, GPUSparseLULinearSolver_Integration){
    rtfem::GPUSparseLULinearSolver<double> gpu_sparse_solver;
    gpu_sparse_solver.PreSolve(*sparse_matrix);
    gpu_sparse_solver.Solve(b.data(), x.data());

    std::cout << "x" << std::endl;
    std::cout << x << std::endl;
}
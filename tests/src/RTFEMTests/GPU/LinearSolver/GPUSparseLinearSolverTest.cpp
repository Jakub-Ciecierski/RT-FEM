#include "RTFEMTests/GPU/LinearSolver/GPUSparseLinearSolverTest.h"

#include <RTFEM/GPU/LinearSolver/GPUSparseLinearSolver.cuh>

void GPUSparseLinearSolverTest::SetUp() {
}

void GPUSparseLinearSolverTest::TearDown() {
}

TEST_F(GPUSparseLinearSolverTest, Test1){
    rtfem::GPUSparseLinearSolver<double> gpu_sparse_solver;
    gpu_sparse_solver.Solve();
}
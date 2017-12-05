#include "RTFEMTests/GPU/LinearSolver/GPULinearSolverTest.h"

#include <RTFEM/GPU/LinearSolver/GPULinearSolver.cuh>

void GPULinearSolverTest::SetUp() {
}

void GPULinearSolverTest::TearDown() {
}

TEST_F(GPULinearSolverTest, Intergration){
    constexpr unsigned int n = 3;
    double A[n*n] = {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0};
    double B[n] = {1.0, 2.0, 3.0};
    double X[n];

    rtfem::GPULinearSolver<double> solver;
    solver.PreSolve(A, n);
    solver.Solve(B, n, X);

    for(unsigned int i = 0; i < n; i++){
        std::cout << "X[" << i << "] = " << X[i] << std::endl;
    }
}
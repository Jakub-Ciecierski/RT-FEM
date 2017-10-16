#include "RTFEMTests/FEM/Solver/FEMSolverTest.h"

#include <RTFEM/FEM/Solver/FEMSolver.h>
#include <RTFEM/FEM/FEMModel.h>
#include <RTFEM/DataTypes.h>
#include <RTFEMTests/Builder/FEMModelBoxBuilder.h>

void FEMSolverTest::SetUp() {
    FEMModelBoxBuilder builder;
    fem_model_ = builder.Create();
}

void FEMSolverTest::TearDown() {
}

TEST_F(FEMSolverTest, FEMSolver_EmptyModel_DisplacemntsAreZero) {
    rtfem::FEMSolver<double> fem_solver;

    auto fem_solver_output = fem_solver.Solve(*fem_model_);

    for(unsigned int i = 0; i < fem_solver_output.displacement.size(); i++){
        EXPECT_EQ(0, fem_solver_output.displacement[i]);
    }
}

TEST_F(FEMSolverTest, FEMSolver_GravityForce_Displacemnts) {
    rtfem::FEMSolver<double> fem_solver;
    fem_model_->SetBodyForce(Eigen::Vector3<double>{
        0, -9.81* 1000000, 0
    });
    auto fem_solver_output = fem_solver.Solve(*fem_model_);

    std::cout << fem_solver_output.displacement << std::endl;
/*
    for(unsigned int i = 0; i < fem_solver_output.displacement.size(); i++){
        EXPECT_EQ(0, fem_solver_output.displacement[i]);
    }*/
}
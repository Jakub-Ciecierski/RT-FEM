#include "RTFEMTests/FEM/BoundaryConditionTest.h"

#include <RTFEM/FEM/Solver/FEMSolver.h>
#include <RTFEM/FEM/Solver/FEMGlobalAssembler.h>
#include <RTFEM/FEM/FEMGeometry.h>
#include <RTFEM/FEM/Vertex.h>
#include "RTFEMTests/Builder/FEMModelSampleBuilder.h"
#include <RTFEM/FEM/BoundaryConditionContainer.h>
#include <RTFEM/FEM/BoundaryCondition.h>

void BoundaryConditionTest::SetUp() {
    fem_model_ = FEMModelSampleBuilder().CreateRandomFEMModel();
}

void BoundaryConditionTest::TearDown() {
}

TEST_F(BoundaryConditionTest, FEMModel_BoundaryAdded_CorrectForceVector) {
    constexpr double boundary_value1 = 1;
    constexpr double boundary_value2 = 2;
    constexpr double boundary_value3 = 3;

    auto size = fem_model_->fem_geometry().vertices.size();
    auto id = fem_model_->fem_geometry().vertices[size - 1]->id();
    fem_model_->boundary_conditions().AddBoundaryCondition(
        rtfem::BoundaryCondition<double>{
            id,
            Eigen::Vector3<double>(boundary_value1,
                                   boundary_value2,
                                   boundary_value3)
        }
    );

    rtfem::FEMGlobalAssembler<double> fem_assembler;

    try{
        auto assembler_data = fem_assembler.Compute(*fem_model_);

        auto force_size = assembler_data.global_force.size();

        EXPECT_EQ(boundary_value3, assembler_data.global_force[force_size - 1]);
        EXPECT_EQ(boundary_value2, assembler_data.global_force[force_size - 2]);
        EXPECT_EQ(boundary_value1, assembler_data.global_force[force_size - 3]);
    }catch(std::exception exception){}

}
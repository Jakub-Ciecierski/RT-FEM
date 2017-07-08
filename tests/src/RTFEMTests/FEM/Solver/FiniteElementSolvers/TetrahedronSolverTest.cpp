#include "RTFEMTests/FEM/Solver/FiniteElementSolvers/TetrahedronSolverTest.h"

#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/Memory/UniquePointer.h>

void TetrahedronSolverTest::SetUp() {
    solver_ = rtfem::make_unique<rtfem::TetrahedronSolver>();

    auto vertex0 = std::make_shared<rtfem::Vertex>(0, rtfem::Vector3(0, 0, 0));
    auto vertex1 = std::make_shared<rtfem::Vertex>(1, rtfem::Vector3(2, 0, 0));
    auto vertex2 = std::make_shared<rtfem::Vertex>(2, rtfem::Vector3(1, 0, 2));
    auto vertex3 = std::make_shared<rtfem::Vertex>(3, rtfem::Vector3(1, 2, 1));
    finite_element_ = std::make_shared<rtfem::TetrahedronFiniteElement>(vertex0,
                                                                        vertex1,
                                                                        vertex2,
                                                                        vertex3);
}

void TetrahedronSolverTest::TearDown() {
}

TEST_F(TetrahedronSolverTest, Solver_GeomtryMatrix_ProperDimensions){
    auto data = solver_->Solve(finite_element_);

    EXPECT_EQ((rtfem::UInt)6, data.geometry_matrix.dimensions().row_count);
    EXPECT_EQ((rtfem::UInt)12, data.geometry_matrix.dimensions().column_count);
}
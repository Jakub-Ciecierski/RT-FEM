#include "RTFEMTests/FEM/Solver/FiniteElementSolvers/TetrahedronSolverTest.h"

#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/Memory/UniquePointer.h>

#include <cmath>

void TetrahedronSolverTest::SetUp() {
    solver_ = rtfem::make_unique<rtfem::TetrahedronSolver>();

    auto vertex0 = std::make_shared<rtfem::Vertex>(0, rtfem::Vector3(1, 0, 2));
    auto vertex1 = std::make_shared<rtfem::Vertex>(1, rtfem::Vector3(1, 2, 1));
    auto vertex2 = std::make_shared<rtfem::Vertex>(2, rtfem::Vector3(0, 0, 0));
    auto vertex3 = std::make_shared<rtfem::Vertex>(3, rtfem::Vector3(2, 0, 0));

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

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest1){
    auto vertex0 = std::make_shared<rtfem::Vertex>(0, rtfem::Vector3(-1.0/2.0, -2.0/9.0, 4.0/11.0));
    auto vertex1 = std::make_shared<rtfem::Vertex>(1, rtfem::Vector3(1.0/16.0, 0, 2.0/13.0));
    auto vertex2 = std::make_shared<rtfem::Vertex>(2, rtfem::Vector3(-1.0/34.0, -1.0/13.0, 1.0/3.0));
    auto vertex3 = std::make_shared<rtfem::Vertex>(3, rtfem::Vector3(-3.0/11.0, -2.0/9.0, -1.0/4.0));
    auto finite_element = std::make_shared<rtfem::TetrahedronFiniteElement>(vertex0,
                                                                            vertex1,
                                                                            vertex2,
                                                                            vertex3);

    rtfem::Matrix expected_jacobian_inverse(4,4);
    expected_jacobian_inverse[0][0] = 0.18591;
    expected_jacobian_inverse[1][0] = 1.54868;
    expected_jacobian_inverse[2][0] = -0.839161;
    expected_jacobian_inverse[3][0] = 0.104569;

    expected_jacobian_inverse[0][1] = -3.65436;
    expected_jacobian_inverse[1][1] = -4.59225;
    expected_jacobian_inverse[2][1] = 7.02344;
    expected_jacobian_inverse[3][1] = 1.22316;

    expected_jacobian_inverse[0][2] = 5.0108;
    expected_jacobian_inverse[1][2] = 14.5185;
    expected_jacobian_inverse[2][2] = -15.3223;
    expected_jacobian_inverse[3][2] = -4.20691;

    expected_jacobian_inverse[0][3] = 0.276165;
    expected_jacobian_inverse[1][3] = -1.70083;
    expected_jacobian_inverse[2][3] = 2.60127;
    expected_jacobian_inverse[3][3] = -1.17661;

    auto jacobian_inverse = solver_->SolveJacobianInverse(finite_element);

    // Round
    for(rtfem::UInt i = 0; i < jacobian_inverse.dimensions().row_count;i++){
        for(rtfem::UInt j = 0; j < jacobian_inverse.dimensions().column_count;j++){
            jacobian_inverse[i][j] = std::floor(jacobian_inverse[i][j] * 100) / 100;
            expected_jacobian_inverse[i][j] = std::floor(expected_jacobian_inverse[i][j] * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest2){
    auto vertex0 = std::make_shared<rtfem::Vertex>(0, rtfem::Vector3(1.0/56.0, -2.0/13.0, -4.0/11.0));
    auto vertex1 = std::make_shared<rtfem::Vertex>(1, rtfem::Vector3(2.0/9.0, 1.0/9.0, -2.0/11.0));
    auto vertex2 = std::make_shared<rtfem::Vertex>(2, rtfem::Vector3(-1.0/9.0, -7.0/15.0, 3.0/7.0));
    auto vertex3 = std::make_shared<rtfem::Vertex>(3, rtfem::Vector3(-1.0/34.0, -4.0/9.0, 2.0/7.0));
    auto finite_element = std::make_shared<rtfem::TetrahedronFiniteElement>(vertex0,
                                                                            vertex1,
                                                                            vertex2,
                                                                            vertex3);

    rtfem::Matrix expected_jacobian_inverse(4,4);
    expected_jacobian_inverse[0][0] = 0.380431;
    expected_jacobian_inverse[1][0] = 0.679071;
    expected_jacobian_inverse[2][0] = 1.95165;
    expected_jacobian_inverse[3][0] = -2.01115;

    expected_jacobian_inverse[0][1] = -3.14555;
    expected_jacobian_inverse[1][1] = 1.23514;
    expected_jacobian_inverse[2][1] = -10.2557;
    expected_jacobian_inverse[3][1] = 12.1661;

    expected_jacobian_inverse[0][2] = -0.102581;
    expected_jacobian_inverse[1][2] = 2.11141;
    expected_jacobian_inverse[2][2] = 6.44379;
    expected_jacobian_inverse[3][2] = -8.45262;

    expected_jacobian_inverse[0][3] = -1.81488;
    expected_jacobian_inverse[1][3] = 1.03481;
    expected_jacobian_inverse[2][3] = 2.13718;
    expected_jacobian_inverse[3][3] = -1.35711;

    auto jacobian_inverse = solver_->SolveJacobianInverse(finite_element);

    // Round
    for(rtfem::UInt i = 0; i < jacobian_inverse.dimensions().row_count;i++){
        for(rtfem::UInt j = 0; j < jacobian_inverse.dimensions().column_count;j++){
            jacobian_inverse[i][j] = std::floor(jacobian_inverse[i][j] * 100) / 100;
            expected_jacobian_inverse[i][j] = std::floor(expected_jacobian_inverse[i][j] * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

TEST_F(TetrahedronSolverTest, Solver_BodyForceAppliedGravity_ProperResult){
    const rtfem::Float g = -9.81;
    rtfem::Vector3 gravity(0, g, 0);
    auto data = solver_->Solve(finite_element_, gravity, rtfem::TractionForce{});

    // Assume volume is calculated correctly.
    auto volume = data.volume;

    rtfem::Matrix expected_node_force_vector(12, 1);
    expected_node_force_vector[0][0] = 0;
    expected_node_force_vector[1][0] = g;
    expected_node_force_vector[2][0] = 0;
    expected_node_force_vector[3][0] = 0;
    expected_node_force_vector[4][0] = g;
    expected_node_force_vector[5][0] = 0;
    expected_node_force_vector[6][0] = 0;
    expected_node_force_vector[7][0] = g;
    expected_node_force_vector[8][0] = 0;
    expected_node_force_vector[9][0] = 0;
    expected_node_force_vector[10][0] = g;
    expected_node_force_vector[11][0] = 0;
    expected_node_force_vector = expected_node_force_vector * (volume/4.0);

    EXPECT_EQ(expected_node_force_vector, data.force_vector);
}

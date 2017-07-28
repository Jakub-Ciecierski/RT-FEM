#include "RTFEMTests/FEM/Solver/FiniteElementSolvers/TetrahedronSolverTest.h"

#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolver.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/Memory/UniquePointer.h>

#include <cmath>

void TetrahedronSolverTest::SetUp() {
    solver_ = rtfem::make_unique<rtfem::TetrahedronSolver<double>>();

    auto vertex0 = std::make_shared<rtfem::Vertex<double>>(0, Eigen::Vector3<double>(1, 0, 2));
    auto vertex1 = std::make_shared<rtfem::Vertex<double>>(1, Eigen::Vector3<double>(1, 2, 1));
    auto vertex2 = std::make_shared<rtfem::Vertex<double>>(2, Eigen::Vector3<double>(0, 0, 0));
    auto vertex3 = std::make_shared<rtfem::Vertex<double>>(3, Eigen::Vector3<double>(2, 0, 0));

    finite_element_ = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertex0,
                                                                                      vertex1,
                                                                                      vertex2,
                                                                                      vertex3);
}

void TetrahedronSolverTest::TearDown() {
}

TEST_F(TetrahedronSolverTest, Solver_GeomtryMatrix_ProperDimensions){
    auto data = solver_->Solve(finite_element_);

    EXPECT_EQ((unsigned int)6, data.geometry_matrix.rows());
    EXPECT_EQ((unsigned int)12, data.geometry_matrix.cols());
}

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest1){
    auto vertex0 = std::make_shared<rtfem::Vertex<double>>(0, Eigen::Vector3<double>(-1.0/2.0, -2.0/9.0, 4.0/11.0));
    auto vertex1 = std::make_shared<rtfem::Vertex<double>>(1, Eigen::Vector3<double>(1.0/16.0, 0, 2.0/13.0));
    auto vertex2 = std::make_shared<rtfem::Vertex<double>>(2, Eigen::Vector3<double>(-1.0/34.0, -1.0/13.0, 1.0/3.0));
    auto vertex3 = std::make_shared<rtfem::Vertex<double>>(3, Eigen::Vector3<double>(-3.0/11.0, -2.0/9.0, -1.0/4.0));
    auto finite_element = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertex0,
                                                                                          vertex1,
                                                                                          vertex2,
                                                                                          vertex3);
    Eigen::Matrix<double, 4, 4> expected_jacobian_inverse;
    expected_jacobian_inverse(0, 0) = 0.18591;
    expected_jacobian_inverse(1, 0) = 1.54868;
    expected_jacobian_inverse(2, 0) = -0.839161;
    expected_jacobian_inverse(3, 0) = 0.104569;

    expected_jacobian_inverse(0, 1) = -3.65436;
    expected_jacobian_inverse(1, 1) = -4.59225;
    expected_jacobian_inverse(2, 1) = 7.02344;
    expected_jacobian_inverse(3, 1) = 1.22316;

    expected_jacobian_inverse(0, 2) = 5.0108;
    expected_jacobian_inverse(1, 2) = 14.5185;
    expected_jacobian_inverse(2, 2) = -15.3223;
    expected_jacobian_inverse(3, 2) = -4.20691;

    expected_jacobian_inverse(0, 3) = 0.276165;
    expected_jacobian_inverse(1, 3) = -1.70083;
    expected_jacobian_inverse(2, 3) = 2.60127;
    expected_jacobian_inverse(3, 3) = -1.17661;

    auto jacobian_inverse = solver_->SolveJacobianInverse(finite_element);

    // Round
    for(unsigned int i = 0; i < jacobian_inverse.rows();i++){
        for(unsigned int j = 0; j < jacobian_inverse.cols();j++){
            jacobian_inverse(i,j) = std::floor(jacobian_inverse(i,j) * 100) / 100;
            expected_jacobian_inverse(i, j) = std::floor(expected_jacobian_inverse(i, j) * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest2){
    auto vertex0 = std::make_shared<rtfem::Vertex<double>>(0, Eigen::Vector3<double>(1.0/56.0, -2.0/13.0, -4.0/11.0));
    auto vertex1 = std::make_shared<rtfem::Vertex<double>>(1, Eigen::Vector3<double>(2.0/9.0, 1.0/9.0, -2.0/11.0));
    auto vertex2 = std::make_shared<rtfem::Vertex<double>>(2, Eigen::Vector3<double>(-1.0/9.0, -7.0/15.0, 3.0/7.0));
    auto vertex3 = std::make_shared<rtfem::Vertex<double>>(3, Eigen::Vector3<double>(-1.0/34.0, -4.0/9.0, 2.0/7.0));
    auto finite_element = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(vertex0,
                                                                                          vertex1,
                                                                                          vertex2,
                                                                                          vertex3);

    Eigen::Matrix<double, 4, 4> expected_jacobian_inverse;
    expected_jacobian_inverse(0, 0) = 0.380431;
    expected_jacobian_inverse(1, 0) = 0.679071;
    expected_jacobian_inverse(2, 0) = 1.95165;
    expected_jacobian_inverse(3, 0) = -2.01115;

    expected_jacobian_inverse(0, 1) = -3.14555;
    expected_jacobian_inverse(1, 1) = 1.23514;
    expected_jacobian_inverse(2, 1) = -10.2557;
    expected_jacobian_inverse(3, 1) = 12.1661;

    expected_jacobian_inverse(0, 2) = -0.102581;
    expected_jacobian_inverse(1, 2) = 2.11141;
    expected_jacobian_inverse(2, 2) = 6.44379;
    expected_jacobian_inverse(3, 2) = -8.45262;

    expected_jacobian_inverse(0, 3) = -1.81488;
    expected_jacobian_inverse(1, 3) = 1.03481;
    expected_jacobian_inverse(2, 3) = 2.13718;
    expected_jacobian_inverse(3, 3) = -1.35711;

    auto jacobian_inverse = solver_->SolveJacobianInverse(finite_element);

    // Round
    for(unsigned int i = 0; i < jacobian_inverse.rows();i++){
        for(unsigned int j = 0; j < jacobian_inverse.cols();j++){
            jacobian_inverse(i, j) = std::floor(jacobian_inverse(i, j) * 100) / 100;
            expected_jacobian_inverse(i, j) = std::floor(expected_jacobian_inverse(i, j) * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

TEST_F(TetrahedronSolverTest, Solver_BodyForceAppliedGravity_ProperResult){
    const double g = -9.81;
    Eigen::Vector3<double> gravity(0, g, 0);
    auto data = solver_->Solve(finite_element_, gravity, rtfem::TractionForce<double>{});

    // Assume volume is calculated correctly.
    auto volume = data.volume;

    Eigen::Matrix<double, 12, 1> expected_node_force_vector;
    expected_node_force_vector(0) = 0;
    expected_node_force_vector(1) = g;
    expected_node_force_vector(2) = 0;
    expected_node_force_vector(3) = 0;
    expected_node_force_vector(4) = g;
    expected_node_force_vector(5) = 0;
    expected_node_force_vector(6) = 0;
    expected_node_force_vector(7) = g;
    expected_node_force_vector(8) = 0;
    expected_node_force_vector(9) = 0;
    expected_node_force_vector(10) = g;
    expected_node_force_vector(11) = 0;
    expected_node_force_vector *= (volume/4.0);

    EXPECT_EQ(expected_node_force_vector, data.force_vector);
}

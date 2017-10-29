#include "RTFEMTests/FEM/Solver/FiniteElementSolvers/TetrahedronSolverTest.h"

#include <RTFEM/FEM/Solver/FiniteElementSolvers/TetrahedronSolver/TetrahedronSolver.h>
#include <RTFEM/FEM/Vertex.h>
#include <RTFEM/FEM/FiniteElements/TetrahedronFiniteElement.h>
#include <RTFEM/Memory/UniquePointer.h>
#include <RTFEM/FEM/Meshing/TriangleMesh.h>

#include <cmath>

void TetrahedronSolverTest::SetUp() {
    solver_ = rtfem::make_unique<rtfem::TetrahedronSolver<double>>();

    vertices_.push_back(std::make_shared<rtfem::Vertex<double>>(0,
                                                                Eigen::Vector3<
                                                                    double>
                                                                    (1, 0, 2)));

    vertices_.push_back(std::make_shared<rtfem::Vertex<double>>(1,
                                                                Eigen::Vector3<
                                                                    double>
                                                                    (1, 2, 1)));
    vertices_.push_back(std::make_shared<rtfem::Vertex<double>>(2,
                                                                Eigen::Vector3<
                                                                    double>(0,
                                                                            0,
                                                                            0)));
    vertices_.push_back(std::make_shared<rtfem::Vertex<double>>(3,
                                                                Eigen::Vector3<
                                                                    double>(2,
                                                                            0,
                                                                            0)));

    triangle_faces_.push_back(rtfem::TriangleFace<double>{1,2,3});
    triangle_faces_.push_back(rtfem::TriangleFace<double>{2,3,0});
    triangle_faces_.push_back(rtfem::TriangleFace<double>{3,0,1});
    triangle_faces_.push_back(rtfem::TriangleFace<double>{0,1,2});

    finite_element_ = std::make_shared<rtfem::TetrahedronFiniteElement<double>>(
        0, 1, 2, 3,
        0, 1, 2, 3);
}

void TetrahedronSolverTest::TearDown() {
}

TEST_F(TetrahedronSolverTest, Solver_GeomtryMatrix_ProperDimensions) {
    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   Eigen::Vector3<double>::Zero());

        EXPECT_EQ((unsigned int) 6, data.geometry_matrix.rows());
        EXPECT_EQ((unsigned int) 12, data.geometry_matrix.cols());
    }catch(std::exception exception){}
}

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest1) {
    auto vertices = std::vector<std::shared_ptr<rtfem::Vertex<double>>>(4);
    vertices[0] = std::make_shared<rtfem::Vertex<double>>(0,
                                                          Eigen::Vector3<double>(
                                                              -1.0 / 2.0,
                                                              -2.0 / 9.0,
                                                              4.0 / 11.0));
    vertices[1] = std::make_shared<rtfem::Vertex<double>>(1,
                                                          Eigen::Vector3<double>(
                                                              1.0 / 16.0,
                                                              0,
                                                              2.0 / 13.0));
    vertices[2] = std::make_shared<rtfem::Vertex<double>>(2,
                                                          Eigen::Vector3<double>(
                                                              -1.0 / 34.0,
                                                              -1.0 / 13.0,
                                                              1.0 / 3.0));
    vertices[3] = std::make_shared<rtfem::Vertex<double>>(3,
                                                          Eigen::Vector3<double>(
                                                              -3.0 / 11.0,
                                                              -2.0 / 9.0,
                                                              -1.0 / 4.0));

    std::vector<rtfem::TriangleFace<double>> triangle_faces;
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 1, 2});
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 1, 3});
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 3, 2});
    triangle_faces.push_back(rtfem::TriangleFace<double>{3, 1, 2});

    auto finite_element =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(
            0, 1, 2, 3,
            0, 1, 2, 3);

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

    auto jacobian_inverse = solver_->SolveJacobianInverse(finite_element,
                                                          vertices);

    // Round
    for (unsigned int i = 0; i < jacobian_inverse.rows(); i++) {
        for (unsigned int j = 0; j < jacobian_inverse.cols(); j++) {
            jacobian_inverse(i, j) =
                std::floor(jacobian_inverse(i, j) * 100) / 100;
            expected_jacobian_inverse(i, j) =
                std::floor(expected_jacobian_inverse(i, j) * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

/**
 * TODO: Should move it to benchmarks test
 */
TEST_F(TetrahedronSolverTest, Solver_JacobianInverse_RandomMathematicaTest2) {
    auto vertices = std::vector<std::shared_ptr<rtfem::Vertex<double>>>(4);

    vertices[0] = std::make_shared<rtfem::Vertex<double>>(0,
                                                          Eigen::Vector3<double>(
                                                              1.0 / 56.0,
                                                              -2.0 / 13.0,
                                                              -4.0 / 11.0));
    vertices[1] = std::make_shared<rtfem::Vertex<double>>(1,
                                                          Eigen::Vector3<double>(
                                                              2.0 / 9.0,
                                                              1.0 / 9.0,
                                                              -2.0 / 11.0));
    vertices[2] = std::make_shared<rtfem::Vertex<double>>(2,
                                                          Eigen::Vector3<double>(
                                                              -1.0 / 9.0,
                                                              -7.0 / 15.0,
                                                              3.0 / 7.0));
    vertices[3] = std::make_shared<rtfem::Vertex<double>>(3,
                                                          Eigen::Vector3<double>(
                                                              -1.0 / 34.0,
                                                              -4.0 / 9.0,
                                                              2.0 / 7.0));

    std::vector<rtfem::TriangleFace<double>> triangle_faces;
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 1, 2});
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 1, 3});
    triangle_faces.push_back(rtfem::TriangleFace<double>{0, 3, 2});
    triangle_faces.push_back(rtfem::TriangleFace<double>{3, 1, 2});

    auto finite_element =
        std::make_shared<rtfem::TetrahedronFiniteElement<double>>(
            0, 1, 2, 3,
            0, 1, 2, 3);

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

    auto jacobian_inverse = solver_->SolveJacobianInverse(
        finite_element,
        vertices);

    // Round
    for (unsigned int i = 0; i < jacobian_inverse.rows(); i++) {
        for (unsigned int j = 0; j < jacobian_inverse.cols(); j++) {
            jacobian_inverse(i, j) =
                std::floor(jacobian_inverse(i, j) * 100) / 100;
            expected_jacobian_inverse(i, j) =
                std::floor(expected_jacobian_inverse(i, j) * 100) / 100;
        }
    }

    EXPECT_EQ(expected_jacobian_inverse, jacobian_inverse);
}

TEST_F(TetrahedronSolverTest, Solver_BodyForceAppliedGravity_ProperResult) {
    const double g = -9.81;
    Eigen::Vector3<double> gravity(0, g, 0);

    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   gravity);

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
        expected_node_force_vector *= (volume / 4.0);

        EXPECT_EQ(expected_node_force_vector, data.force_vector);
    }catch(std::exception exception){}
}

TEST_F(TetrahedronSolverTest, Solver_TractionForceAppliedToFace1_ProperResult) {
    unsigned int triangle_face_index = 0;
    triangle_faces_[triangle_face_index].traction_force = -9.81;

    auto index1 = triangle_face_index*3 + 0;
    auto index2 = triangle_face_index*3 + 0;
    auto index3 = triangle_face_index*3 + 0;

    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   Eigen::Vector3<double>::Zero());
        for(unsigned int i = 0; i < 12; i++){
            if(i == index1 || i == index2 || i == index3){
                EXPECT_EQ(0, data.force_vector[i]);
            }else{
                EXPECT_NE(0, data.force_vector[i]);
            }
        }

    }catch(std::exception exception){}
}

TEST_F(TetrahedronSolverTest, Solver_TractionForceAppliedToFace2_ProperResult) {
    unsigned int triangle_face_index = 1;
    triangle_faces_[triangle_face_index].traction_force = -9.81;

    auto index1 = triangle_face_index*3 + 0;
    auto index2 = triangle_face_index*3 + 0;
    auto index3 = triangle_face_index*3 + 0;

    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   Eigen::Vector3<double>::Zero());
        for(unsigned int i = 0; i < 12; i++){
            if(i == index1 || i == index2 || i == index3){
                EXPECT_EQ(0, data.force_vector[i]);
            }else{
                EXPECT_NE(0, data.force_vector[i]);
            }
        }

    }catch(std::exception exception){}
}

TEST_F(TetrahedronSolverTest, Solver_TractionForceAppliedToFace3_ProperResult) {
    unsigned int triangle_face_index = 2;
    triangle_faces_[triangle_face_index].traction_force = -9.81;

    auto index1 = triangle_face_index*3 + 0;
    auto index2 = triangle_face_index*3 + 0;
    auto index3 = triangle_face_index*3 + 0;

    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   Eigen::Vector3<double>::Zero());
        for(unsigned int i = 0; i < 12; i++){
            if(i == index1 || i == index2 || i == index3){
                EXPECT_EQ(0, data.force_vector[i]);
            }else{
                EXPECT_NE(0, data.force_vector[i]);
            }
        }

    }catch(std::exception exception){}
}

TEST_F(TetrahedronSolverTest, Solver_TractionForceAppliedToFace4_ProperResult) {
    unsigned int triangle_face_index = 3;
    triangle_faces_[triangle_face_index].traction_force = -9.81;

    auto index1 = triangle_face_index*3 + 0;
    auto index2 = triangle_face_index*3 + 0;
    auto index3 = triangle_face_index*3 + 0;

    try{
        auto data = solver_->Solve(finite_element_,
                                   vertices_,
                                   triangle_faces_,
                                   Eigen::Vector3<double>::Zero());
        for(unsigned int i = 0; i < 12; i++){
            if(i == index1 || i == index2 || i == index3){
                EXPECT_EQ(0, data.force_vector[i]);
            }else{
                EXPECT_NE(0, data.force_vector[i]);
            }
        }

    }catch(std::exception exception){}
}
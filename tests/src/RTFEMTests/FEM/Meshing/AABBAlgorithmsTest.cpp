#include "RTFEMTests/FEM/Meshing/AABBAlgorithmsTest.h"

#include <RTFEM/FEM/Meshing/TriangleMesh.h>
#include "RTFEM/FEM/Meshing/AABBAlgorithms.h"
#include "RTFEM/FEM/Meshing/AABB.h"
#include <RTFEM/Memory/UniquePointer.h>

void AABBAlgorithmsTest::SetUp() {
    triangle_mesh_ = rtfem::make_unique<rtfem::TriangleMesh<float>>();
    triangle_mesh_->triangles = {
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{-2.0f, 9.0f, 0.0f},
                               Eigen::Vector3<float>{-3.0f, 1.0f, 0.0f},
                               Eigen::Vector3<float>{5.0f, 2.0f, 3.0f}},
        rtfem::TriangleFaceWithPoints<float>{Eigen::Vector3<float>{-2.0f, 0.5f, 2.0f},
                               Eigen::Vector3<float>{-3.0f, -1.0f, 1.0f},
                               Eigen::Vector3<float>{6.0f, 2.0f, 2.0f}},
    };

}

void AABBAlgorithmsTest::TearDown() {
}

TEST_F(AABBAlgorithmsTest, FindMin_ProperXMin) {
    constexpr float expected_min = -3.0f;
    EXPECT_EQ(expected_min,
              rtfem::FindMin(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::X));
}

TEST_F(AABBAlgorithmsTest, FindMin_ProperYMin) {
    constexpr float expected_min = 1.0f;
    EXPECT_EQ(expected_min,
              rtfem::FindMin(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::Y));
}

TEST_F(AABBAlgorithmsTest, FindMin_ProperZMin) {
    constexpr float expected_min = 0.0f;
    EXPECT_EQ(expected_min,
              rtfem::FindMin(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::Z));
}

TEST_F(AABBAlgorithmsTest, FindMax_ProperXMax) {
    constexpr float expected_max = 5.0f;
    EXPECT_EQ(expected_max,
              rtfem::FindMax(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::X));
}

TEST_F(AABBAlgorithmsTest, FindMax_ProperYMax) {
    constexpr float expected_max = 9.0f;
    EXPECT_EQ(expected_max,
              rtfem::FindMax(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::Y));
}

TEST_F(AABBAlgorithmsTest, FindMax_ProperZMax) {
    constexpr float expected_max = 3.0f;
    EXPECT_EQ(expected_max,
              rtfem::FindMax(triangle_mesh_->triangles[0],
                             rtfem::AABBCoordinate::Z));
}

TEST_F(AABBAlgorithmsTest, CreateAABB) {
    auto aabb = rtfem::CreateAABB<float>(*triangle_mesh_);

    EXPECT_EQ(6, aabb.Max()[0]);
    EXPECT_EQ(9, aabb.Max()[1]);
    EXPECT_EQ(3, aabb.Max()[2]);

    EXPECT_EQ(-3, aabb.Min()[0]);
    EXPECT_EQ(-1, aabb.Min()[1]);
    EXPECT_EQ(0, aabb.Min()[2]);
}